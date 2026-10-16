# Lab Test ITS2 PCR-cycle analysis

연구실에서 수행한 human skin ITS2 시험 데이터의 전처리 조건, 음성대조군 오염 양상, PCR cycle 효과 및 taxonomy 안정성을 점검하는 프로젝트이다.

> 연구 진행 기록과 상세 메모: [Notion — Test ITS2 PCR cycle](https://app.notion.com/p/wait-a-minute/Test-ITS2-PCR-cycle-3fa27c45489880daa08ccb8613bd9fb8?source=copy_link)

## 연구 질문

이 프로젝트는 다음 질문에서 시작되었다.

1. 현재 연구실 workflow로 피부 진균 ITS2 데이터를 안정적으로 얻을 수 있는가?
2. True skin sample과 negative control의 alpha diversity, beta diversity 및 taxonomic composition은 어떻게 다른가?
3. Air swab, clean-bench swab, empty-tube control에는 어떤 오염 profile이 존재하는가?
4. PCR 조건과 DNA 농도에 따라 검출 ASV와 잠재적 contaminant가 어떻게 변하는가?
5. Cutadapt와 ITSxpress, paired-end와 single-end, DADA2 `maxEE` 설정이 결과에 미치는 영향은 무엇인가?
6. Negative control의 read 수가 적지 않은 이유가 실제 환경·시약 오염인지, PCR 증폭인지, trimming 및 분류 과정의 산물인지 구분할 수 있는가?

## 데이터와 비교 조건

- True sample: human skin sample
  - 대상자와 skin site를 구분하여 해석
  - 동일 biological sample에 여러 PCR condition이 존재
- Controls
  - Air swab
  - Clean-bench swab
  - Empty-tube control
- 주요 preprocessing 비교
  - QC20 및 QC30
  - DADA2 `maxEE 2/2`, `4/4`, `8/8`
  - Cutadapt primer trimming
  - ITSxpress ITS2 extraction
  - Paired-end 및 forward single-end

Rarefaction curve는 library depth의 충분성과 다양성 포화 추이를 확인하는 용도이다. 임의의 rarefaction cutoff로 sample을 사전 제거하지 않으며, taxonomic composition은 rarefied table이 아닌 raw count 또는 total-sum-scaled relative abundance를 사용한다.

## 지금까지의 분석 과정과 고민

### 1. Cutadapt와 ITSxpress의 결과가 예상보다 크게 달랐다

처음에는 ITSxpress가 보존된 rRNA flank를 제거하고 ITS2만 남기므로 taxonomy resolution을 높일 것으로 기대했다. 실제로 일부 Cutadapt ASV 안에 더 짧은 ITSxpress ASV가 정확히 포함되는 `exact containment` 관계가 확인되었다.

그러나 read mode와 입력 구성에 따라 ITSxpress의 결과가 크게 달라졌다.

- 초기 single/paired 비교에서는 Aureobasidium-like family가 Cutadapt single, Cutadapt paired 및 ITSxpress single에서 회수되었지만 ITSxpress paired에서 주로 소실되었다.
- complete-forward 비교에서는 Cutadapt가 약 1,186,913 reads와 1,855 ASVs를 보존한 반면 ITSxpress는 약 115,907 reads와 211 ASVs만 보존했다.
- complete-forward sample별 ITSxpress/Cutadapt read recovery 중앙값은 약 6%였다.
- 이 과정에서 식물 유래로 확인된 ASV뿐 아니라 `Malassezia restricta`, `M. globosa`, `Aureobasidium` 계열 신호도 함께 크게 줄었다.

따라서 ITSxpress 통과 여부를 곧바로 “진짜 진균 여부”로 해석할 수 없다. 현재 증거는 incomplete 또는 특정 구조의 forward read에서 ITSxpress HMM boundary 인식이 taxon- 및 sequence-dependent selection을 일으킬 가능성을 시사한다.

### 2. Fungi_k와 Ascomycota_p는 단순한 짧은 ASV 문제가 아니었다

Cutadapt에서 `Fungi_k` 또는 `Ascomycota_p`까지만 분류된 ASV와 ITSxpress ASV를 비교한 결과, ITSxpress 서열이 더 긴 Cutadapt ASV의 내부에 완전히 동일하게 포함되는 사례가 확인되었다.

이때 Cutadapt에만 남아 있는 외부 구간을 처음에는 primer 잔존 또는 반복서열로 의심했다. 하지만 서열 구조만으로는 다음을 구분할 수 없다.

- 실제 5.8S/28S flank
- 제거되지 않은 primer 또는 primer-adjacent sequence
- low-complexity/repetitive sequence
- 비표적 진핵생물 amplicon

따라서 raw read boundary, primer 위치, reference alignment 및 BLAST를 함께 확인해야 한다.

### 3. `Aureobasidium tremulum` 표기는 보수적으로 수정해야 한다

Cutadapt에서 `Aureobasidium tremulum`으로 분류된 ASV가 특히 한 대상자의 True sample에서 많이 관찰되었고, 농도와 양의 상관관계를 보이는 탐색적 결과가 있었다.

대표 full ASV와 ITS2 core의 BLAST는 모두 `Aureobasidium`을 지지했다. 하지만 여러 종이 동일한 query coverage와 identity로 동률이었기 때문에 현재 근거로 `A. tremulum` 종 수준을 확정할 수 없다.

현재 권장 표현은 다음과 같다.

- `Aureobasidium-like ASV`
- `Aureobasidium spp.`
- `Aureobasidium-associated sequence family`

특정 대상자에서의 증가가 실제 피부 mycobiome 차이인지, DNA 농도/PCR 증폭 효과인지, 또는 trimming·분류 편향인지는 아직 해결되지 않았다.

### 4. KSY sample에서는 Malassezia 신호가 `Fungi_k`로 남았을 가능성을 의심했다

KSY 피부 sample에서는 예상했던 Malassezia가 genus/species 수준으로 충분히 분류되지 않는 대신 고풍부 `Fungi_k` ASV가 관찰되었다. 따라서 일부 Malassezia 서열이 낮은 taxonomic rank에서 멈춰 `Fungi_k`로 표시되었을 가능성을 먼저 의심했다.

다만 실제 서열 확인 결과, 대표적인 281-bp `Fungi_k` ASV `5c8f111a2f4a8e9e41a1cb58ebc985f6`는 BLAST에서 281/281 bp, 100% identity로 `Raphanus sativus`를 지지했다. 따라서 적어도 이 ASV는 미분류 Malassezia가 아니라 식물 유래 off-target sequence이다.

이 결과가 KSY의 모든 `Fungi_k` ASV가 Raphanus라는 뜻은 아니다. 다른 `Fungi_k` ASV 중에는 reference 또는 서열 구간의 한계로 Malassezia가 낮은 rank에 머문 경우가 있을 수 있으므로, KSY에서 풍부한 `Fungi_k` ASV를 개별적으로 BLAST하고 broad-eukaryote classifier 결과와 대조해야 한다.

### 5. maxEE 2/2와 8/8은 핵심 차이를 설명하지 못했다

ITSxpress paired 결과에서 `maxEE 2/2`와 `8/8`을 비교했을 때 전체 community structure와 주요 species profile은 거의 같았다. `maxEE 8/8`에서만 나타나는 ASV는 대부분 극저빈도·단일 sample ASV였고 새로운 주요 species signal을 만들지 않았다.

따라서 현재 관찰되는 Cutadapt–ITSxpress 차이와 주요 taxon 소실은 `maxEE` 설정만으로 설명하기 어렵다.

### 6. 분석 단위와 시각화 원칙을 정리했다

- Skin sample은 대상자와 site를 함께 고려한다.
- PCR condition은 지정된 순서의 factor로 유지한다.
- Alpha diversity는 Observed ASV, Shannon 및 read count를 서로 다른 plot으로 제시한다.
- Alpha diversity의 대상자별 panel은 서로 다른 row에 배치한다.
- Beta diversity는 Aitchison distance와 PCoA를 중심으로 보고 paired method effect를 함께 평가한다.
- Taxonomic composition은 relative abundance와 raw read count를 모두 확인한다.
- Heatmap annotation은 필요한 변수만 사용하고 중복 정보를 늘리지 않는다.
- ASV length `<=100 bp` 분석은 누적 abundance가 100 이하인 ASV와 구분한다.
- Focal taxon 및 control-source ternary plot은 겹치지 않도록 taxon/조건별로 분리한다.

## 현재까지 비교적 확실한 결론

1. Preprocessing 방법은 read retention, ASV 수 및 taxonomic label에 큰 영향을 준다.
2. `maxEE 2/2`와 `8/8` 차이는 Cutadapt–ITSxpress 차이에 비해 작다.
3. ITSxpress paired 결과에서 특정 fungal signal의 소실을 biological absence로 해석하면 안 된다.
4. complete-forward에서도 ITSxpress가 매우 많은 read를 제거했으므로 ITSxpress 결과만을 primary community profile로 사용하는 것은 위험하다.
5. Cutadapt는 진균 read를 더 많이 보존하지만 식물 등 non-fungal off-target도 포함할 수 있다.
6. `Raphanus sativus` ASV의 존재는 non-fungal sequence를 별도로 분류하고 제거해야 함을 보여준다.
7. classifier의 species label만으로 중요한 ASV의 종을 확정하지 말고, sequence family와 BLAST 근거를 함께 확인해야 한다.
8. Control에서 read가 검출되는 사실만으로 contamination source나 실험 실패를 확정할 수 없다. Control type, 농도, PCR condition 및 ASV identity를 함께 봐야 한다.

## 아직 풀리지 않은 문제

### Biological question

- 특정 대상자의 True sample에서 Aureobasidium-like ASV가 증가한 이유는 무엇인가?
- 관찰된 농도 상관성이 biomass 효과인지 PCR amplification 효과인지 구분 가능한가?
- Control과 True에서 검출되는 Malassezia/Cutaneotrichosporon 계열이 동일 ASV인가, 또는 근연하지만 다른 ASV인가?
- Air swab, clean-bench swab 및 empty-tube control 사이에 반복적으로 공유되는 core contaminant가 존재하는가?
- Negative control의 높은 read 수는 환경 채취, extraction, PCR, index hopping/cross-talk 중 어느 단계와 더 일치하는가?

### Bioinformatic question

- 초기 ITSxpress single에서는 회수된 Aureobasidium-like signal이 complete-forward ITSxpress에서는 다시 크게 감소한 이유는 무엇인가?
- ITSxpress가 인식하지 못하는 ASV에는 공통적인 길이, flank, GC, low-complexity 또는 primer-boundary pattern이 있는가?
- Cutadapt-only Malassezia와 Aureobasidium ASV는 BLAST와 broad-eukaryote taxonomy에서도 진균으로 지지되는가?
- UNITE classifier의 species label과 NCBI nt BLAST가 불일치하는 ASV를 어떤 보수적 규칙으로 처리할 것인가?
- Paired-end merge/trim 실패와 single-end 정보 손실 중 어느 쪽이 community distortion에 더 큰 영향을 주는가?

## 다음 분석 계획

### 우선순위 1: broad-eukaryote taxonomy 재분류

Cutadapt complete-forward DADA2 ASV를 primary 대상으로 유지하고, UNITE+INSD All Eukaryotes classifier로 전체 ASV를 재분류한다.

권장 pretrained classifier:

```text
UNITE v10.0 2025-02-19
dynamic clustering
all eukaryotes
scikit-learn 1.4.2
unite_ver10_dynamic_all_19.02.2025-Q2-2024.10.qza
```

QIIME 2 2024.10, 2025.4, 2025.7 및 2025.10 계열에서 호환성을 확인한 뒤 사용한다. 식물·동물 등 확실한 non-fungal ASV를 별도 flag하고, 곧바로 삭제하기 전에 abundance와 sample distribution을 audit한다.

### 우선순위 2: high-impact ASV batch BLAST

다음 ASV를 우선 검증한다.

- Cutadapt-only 고풍부 ASV
- `Fungi_k`, `Ascomycota_p` 및 species 미분류 ASV
- True 또는 특정 대상자에서 우점하는 ASV
- Control에서 반복 검출되는 ASV
- Cutadapt와 ITSxpress의 차이를 크게 만드는 ASV

평가 항목은 query coverage, percent identity, alignment length, top-hit 간 동률 여부 및 fungal/non-fungal 일관성이다. BLAST top hit 하나만으로 species를 확정하지 않는다.

### 우선순위 3: ITSxpress loss mechanism 분석

- ASV length와 ITSxpress recovery의 관계
- GC content와 recovery의 관계
- left/right flank 길이와 recovery의 관계
- low-complexity 및 homopolymer pattern
- primer/partial-primer match
- taxon별 recovery fraction
- raw read에서 실제 boundary 위치

### 우선순위 4: Control–True 및 contamination 분석

- Control type별 prevalence와 raw abundance
- True 대비 Control enrichment
- DNA concentration을 이용한 decontam frequency method
- Control label을 이용한 prevalence method
- 두 방법 및 threshold에 따른 contaminant call sensitivity
- 동일 PCR option 안에서 세 control source의 ASV sharing과 ternary profile
- contaminant 제거 전후 alpha, Aitchison beta 및 taxonomic composition 비교

반복 PCR 조건과 대상자 구조를 무시한 단순 독립검정은 피하고, 결과는 exploratory/sensitivity analysis로 구분한다.

### 우선순위 5: 최종 primary/sensitivity dataset 결정

현재 잠정안은 다음과 같다.

- Primary: Cutadapt primer-trimmed complete-forward → DADA2
- Taxonomy cleanup: UNITE+INSD All Eukaryotes + 주요 ASV BLAST
- Sensitivity: ITSxpress single/paired 결과
- Taxonomic composition: raw count 및 TSS relative abundance
- Alpha/Beta: 명시적인 library-depth 검토 후 수행
- Contamination: 제거 전후 결과를 함께 제시

## 저장소 구성

- `01_*.Rmd`: preprocessing, method comparison, sequence/BLAST 검증
- `02_*.Rmd`: basic ecological analysis 및 Control–True 비교
- `Results/**.R`: 재현 가능한 진단 분석 코드
- `benchmark_results/**/versions*.yml`: workflow software version 기록
- `Lab_Test_ITS2_2061015.Rproj`: RStudio project 설정

## 데이터 공개 범위

이 저장소에는 분석 코드와 환경 기록만 포함한다. Raw FASTQ, QIIME 2 artifact, phyloseq object, sample-level metadata, 렌더된 HTML 및 분석 산출물은 포함하지 않는다. 연구가 진행 중이고 human-derived sample 정보를 포함하므로 현재 저장소는 private 운영을 권장한다.

공개 전에는 절대경로, 대상자 식별 가능 정보, 중간 해석 및 불필요한 초안을 다시 점검하고, 익명화된 예제 metadata와 재현 방법을 별도로 제공한다.
