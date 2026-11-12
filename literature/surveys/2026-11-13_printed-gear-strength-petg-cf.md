# FDM 印刷歯車の歯元曲げ強度: PETG・PETG-CF/GF・PLA-CF/GF の許容応力を決めるための文献調査

- 作成日: 2026-11-13
- survey_slug: `printed-gear-strength-petg-cf` (依頼主指定のファイル名に合わせた)
- 深度: focused (採録 27 本、うち歯車そのものの試験 9 本)
- 実行形態: 自動実行 (対話確認なし)。Map と Hub 選定は確認を取らずに進めた (Methodology 参照)

## 0. この文書の目的

小型の段付き遊星減速機 (モジュール 0.9 と 1.1、歯数 12〜54、平歯とやまば歯、0.4mm ノズルの FDM 印刷) は、歯元曲げ応力の許容値を仮に 10 MPa と置いて歯幅を決めている (README.md 191〜192 行)。応力は修正 Lewis 式 σ = Ft/(b·π·y·m)、y = 0.154 − 0.912/z、速度係数 Kv = 1、入力ピーク 0.4 N·m で計算する。この許容値には根拠がない。許容値は歯幅、つまり減速機の高さを決める。

本調査は、PETG、PETG-CF (炭素繊維入り)、PETG-GF (ガラス繊維入り)、PLA-CF、PLA-GF について、印刷歯車の強度データと設計規則を集め、上の Lewis 式と同じ基準で使える許容応力を導く。素の PLA とナイロンは比較対象としてのみ扱う。

### 研究課題 (RQ)

- RQ1: 印刷歯車の歯の静的破壊強度と疲労強度 (S–N) はどれだけ測られているか。条件 (材料、印刷設定、モジュール、歯幅、試験機、応力算出法) は何か。
- RQ2: 印刷歯車の許容曲げ応力は、どの安全率と繰返し数から導かれているか。
- RQ3: 繊維充填、印刷条件 (向き、外周数、充填率、層高さ、温度)、焼きなましは歯の強度をどう変えるか。
- RQ4: VDI 2736 や ISO 6336 を FDM 歯車に当てた例と補正係数、摩耗と温度の限界、印刷遊星減速機の実測破損トルクはどうか。

### 採録基準

- 採録: 査読論文と主要な会議論文・プレプリント。FDM/FFF 印刷歯車の強度・寿命試験、対象材料の印刷試験片の曲げ・疲労・焼きなまし試験、印刷減速機の実機報告。
- 除外: 販売元のブログ、出典のない寿命数値、SLS (粉末焼結) のみの研究 (比較で言及する場合は明記)。

## 1. 要点 (先に結論)

1. **PETG・PETG-CF/GF・PLA-CF/GF の印刷歯車について、歯元応力と寿命の組 (S–N) を測った論文は見つからなかった。** 最も近いのは、FFF 印刷の PLA・ABS・PA66・PA66 連続炭素繊維の歯車を VDI 2736 の方法で疲労試験した Dennig ら 2020 である。
2. その試験で、2 外周・充填 22% の ABS 歯車は、VDI の歯元応力 σF ≈ 10〜20 MPa で 10⁴〜10⁵ 回、≈ 6〜9 MPa でも 3〜5×10⁵ 回で歯元から折れた (Dennig2020 p11 図 13、p13)。PA66 は σF = 12 MPa で 2×10⁶ 回まで壊れなかった (同 p14 図 16)。PLA は 1500 rpm で歯が溶けた (同 p12)。
3. PETG の印刷試験片 (45° ラスタ、充填 100%) の引張–引張疲労では、最大応力 10.2 MPa で 0.8〜1.9×10⁵ 回で破断した (Martins2024 表 2)。同じ表から本調査で当てはめると、10⁵ 回で ≈ 11 MPa、10⁶ 回で ≈ 6 MPa になる。
4. 繊維入り PETG の強度上昇は試験片で +8〜24% と小さい (Valvez2022 表 4・5、Rihar2026 要旨)。焼きなましで強度が下がった例もある (Economides2024 §3)。繊維入り材のカタログ強度が高くても、印刷歯車では先に割れた例がある (Yoshida2023 p4–5)。
5. 以上から、Lewis 式・Kv = 1 の基準での推奨許容値は、**PETG 6 MPa、PETG-CF/GF 6 MPa、PLA-CF/GF 6 MPa** (いずれもピーク荷重 10⁵ 回想定、安全率 1.5、確度は低〜中) とした。現在の 10 MPa は、PETG が 10⁵ 回前後で壊れる応力とほぼ同じで、安全率がほぼ 1 である。導出は §7。

## 2. Survey Findings

### thesis (この分野の根本的な未解決問題)

印刷歯車の強度は「材料」ではなく「材料 × 印刷条件 × 歯の形 × 温度」の組で決まる。しかし規格 (VDI 2736) の材料データは射出成形の PA66 と POM に限られ、印刷材のデータはない (Dennig2020 p4)。射出成形 PA6 でさえ、実測の許容歯元応力は VDI 2736 の PA66 値の 38〜72% だった (Pogacnik2015 p10 表 7)。印刷材では、試験片の強度から歯の強度を予測する換算が確立していない。そのため、許容値は個別の歯車試験で決めるしかない、というのが現状である。

### foundation (分野が依拠する技術要素)

- **応力の算出**: VDI 2736 の歯元応力 σF = YFa·YSa·Yε·Ft/(b·mn) (Dennig2020 p4 式 1)。歯形係数 YFa、応力集中係数 YSa、かみ合い率係数 Yε を掛ける。Lewis 式は Yε も応力集中も含まず、歯先 1 点に全荷重がかかると仮定する。
- **試験法**: 鋼歯車と組ませた乾式の運転試験で、各荷重 3 個以上、10⁶〜2×10⁶ 回で打ち切る (Dennig2020 p11)。荷重を 2×10⁵ 回ごとに 0.05 N·m 上げる段階試験で寿命試験を短縮する方法もある (Pogacnik2015 §3.3)。
- **故障モード**: 歯元折損、溶融、摩耗、ピッチング、層間はく離、変形 (Dennig2020 p3 図 2)。印刷歯車ではこれに層の境界での割れが加わる (Yoshida2023 p4–5、Hriberšek2024 §3)。

### progress (解決されてきたことの流れ)

- 2015: 射出成形 PA6 の寿命試験で、VDI 2736 の許容値が楽観的すぎる可能性が示された (Pogacnik2015 p10)。
- 2020: FFF 歯車を VDI 2736 の試験法で初めて系統的に試験した。PLA は溶融、ABS は歯元折損で壊れ、PA66 系は長寿命だった (Dennig2020)。
- 2022: PETG・PLA・ABS の印刷歯車を鋼歯車と組ませた寿命試験で、PETG が最も摩耗が少なく、最も長持ちした (Tunalioglu2022 §3.1–3.2)。
- 2023: 印刷遊星減速機 (Onyx・POTICON) で出力 19.6 N·m を 46 時間出した。一方、カタログ強度の高い材料の太陽歯車が先に割れた (Yoshida2023)。
- 2024: 印刷 CF-PA (Onyx) 歯車の寿命を 3 荷重で測った (Hriberšek2024。のちに「懸念表明」が付いた)。PETG と PLA の印刷試験片の S–N 曲線が出た (Martins2024)。
- 2025〜2026: PETG・PLA・ASA・ABS 歯車の静的破壊比較 (Nguyen2025)。PETG と PETG+CF、PLA と PLA+CF の充填条件の比較 (Rihar2026、Hozdić2025)。PETG の焼きなまし条件 (Spina2026、Zisopol2026)。

### gap (構造的に未解決の問題)

- G1: **PETG・PETG-CF/GF・PLA-CF/GF の歯車の S–N データがない。** 本調査の範囲では 1 件も見つからなかった。
- G2: **小モジュール (m ≤ 1.5) の印刷歯車の疲労データがほぼない。** 歯車疲労の唯一の系統試験は m = 2 である (Dennig2020 p5 表 1)。Hriberšek2024 は m = 1 だが懸念表明付きである。
- G3: **低速・間欠負荷の条件がない。** 既存の歯車試験は 900〜1500 rpm の連続運転で、温度上昇が故障を左右している (Dennig2020 p12、Tunalioglu2022 §2.2、Hriberšek2024 §3)。遠隔操作やロボット関節のような低速でピークが間欠的に入る使い方とは条件が違う。
- G4: **試験片強度と歯の強度の換算がない。** 繊維入り材は試験片の強度が上がっても歯では上がらないことがある (Yoshida2023 p4–5)。
- G5: **やまば歯の印刷歯車の疲労データがない。** 静的破壊の比較はある (Nguyen2025 §3.3.1)。

## 3. Concept Matrix

列の意味: 歯車 = 歯車そのものの試験、疲労 = 繰返し数を測った、静的 = 静的破壊、対象材 = PETG/CF/GF/PLA-CF/GF を含む、印刷条件 = 向き・充填・外周などを変えた、規格 = VDI 2736/ISO 6336/Lewis を使った、温度 = 歯の温度を測った・論じた、実機 = 減速機・アクチュエータ。● = 該当、○ = 部分的、空欄 = 該当なし。

| 論文 | 歯車 | 疲労 | 静的 | 対象材 | 印刷条件 | 規格 | 温度 | 実機 |
|---|---|---|---|---|---|---|---|---|
| Dennig2020 (FFF 歯車の実現性) | ● | ● |  |  | ○ (外周 2→4) | ● VDI | ● |  |
| Pogacnik2015 (射出 PA6 加速試験) | ● | ● |  |  |  | ● VDI | ● |  |
| Hriberšek2024 (印刷 CF-PA 歯車) | ● | ● |  | ○ (CF) | ○ (流量) | ○ VDI 摩耗 | ● |  |
| Tunalioglu2022 (PETG 等歯車寿命) | ● | ○ (相対寿命) |  | ● PETG |  |  |  |  |
| Nguyen2025 (印刷歯車の破壊試験) | ● |  | ● | ● PETG | ● 充填・モジュール |  |  |  |
| Kalani2022 (PA12 歯車の熱と摩耗) | ● | ○ |  |  |  |  | ● |  |
| Zhang2020 (印刷歯車のパラメータ) | ● | ○ (摩耗) |  |  | ● |  |  |  |
| Ciobanu2025 (FDM 歯車の摩耗) | ● | ○ |  |  | ● 充填 |  |  |  |
| Pisula2021 (FDM 歯車の摩耗) | ● | ○ |  |  |  |  |  |  |
| Zisopol2022 (PLA 平歯車の剛性) | ● |  | ○ |  | ● |  |  |  |
| Dennig2019 (PEEK 連続 CF 歯車) | ● | ● |  |  |  | ● VDI | ● |  |
| Pandian2022 (SLS PA12 歯の向き) | ● | ● |  |  | ● 向き |  |  |  |
| Martins2024 (PLA・PETG の S–N) |  | ● |  | ● PETG | ○ 向き |  |  |  |
| Dolzyk2019 (PETG ラスタ角疲労) |  | ● |  | ● PETG | ● ラスタ角 |  |  |  |
| Kargar2023 (CF-PLA 疲労) |  | ● |  | ● PLA-CF | ● 充填・ラスタ角 |  |  |  |
| CF-PETG S–N (Springer 2023) |  | ● |  | ● PETG-CF |  |  |  |  |
| Valvez2022 (PETG/CF/KF 最適化) |  |  | ● 曲げ | ● PETG-CF | ● |  |  |  |
| Economides2024 (CF-PETG と焼きなまし) |  |  | ● | ● PETG-CF | ● 充填・焼きなまし |  |  |  |
| Rihar2026 (PETG と PETG+CF) |  |  | ● | ● PETG-CF | ● 充填 |  |  |  |
| Hozdić2025 (PLA と PLA+CF) |  |  | ● | ● PLA-CF | ● 充填 |  |  |  |
| Spina2026 (PET/PETG の焼きなまし) |  |  | ● | ● PETG | ● 焼きなまし |  |  |  |
| Zisopol2026 (PETG 焼きなまし曲げ) |  |  | ● 曲げ | ● PETG | ● 層高さ・充填 |  |  |  |
| Yoshida2023 (印刷遊星 QDD) | ● | ○ (連続運転) |  | ○ (繊維入り) | ○ |  |  | ● |
| Chi2025 (Berkeley Humanoid Lite) | ○ (サイクロイド) | ○ (60 h) |  |  |  |  |  | ● |
| Jenks2025 (OpenSEA) | ○ |  |  |  |  |  |  | ● |
| Bezzini2026 (複合遊星と印刷軸受) | ● |  |  |  |  |  |  | ● |
| Singh2025 (外歯・内歯遊星の比較) |  |  |  |  |  | ● Lewis |  | ○ |

疎な箇所: 「歯車 × 疲労 × 対象材」の 3 列が同時に埋まる行は 1 本もない (Tunalioglu2022 は相対寿命のみで応力値がない)。これが G1 である。

## 4. Quantitative Trends

- 年別: 2015 年 1 本、2019〜2020 年 4 本、2021〜2023 年 7 本、2024 年 5 本、2025〜2026 年 10 本。印刷材の試験片研究は 2024 年以降に増えたが、歯車の疲労試験は増えていない。
- 掲載先: MDPI Polymers が 9 本で最多。歯車試験の専門会議 (Dennig2019・2020) と IEEE/SICE SII (Yoshida2023) が少数ある。
- 試験の種別: 歯車の運転試験 9 本、試験片試験 10 本、実機報告 5 本、設計手法 1 本。

## 5. Hub Papers (深読み)

| # | 論文 | 年 | 掲載先 | 選んだ理由 |
|---|---|---|---|---|
| 1 | Dennig, Zumofen, Kirchheim「Feasibility investigation of gears manufactured by fused filament fabrication」 | 2020 | AMPA 2020 (Springer) | B+C: 「印刷条件」と「VDI 規格」の 2 クラスタを 1 本で結ぶ唯一の FFF 歯車疲労試験。thesis と推奨値の主たる根拠 |
| 2 | Martins ら「Mechanical Properties of Additively Manufactured Polymeric Materials—PLA and PETG」 | 2024 | Polymers | C: PETG と PLA の S–N を表で出した唯一の資料。推奨値の材料側の根拠 |
| 3 | Pogacnik, Tavčar「Accelerated testing method for polymer gears」 | 2015 | 会議論文 (掲載先未確認) | C: 規格値と実測許容値の比 (38〜72%) を示し、規格値をそのまま使えない根拠になる |
| 4 | Yoshida ら「Experimental Evaluation of a QDD Actuator with a 3D-printed Planetary Gear Reducer」 | 2023 | IEEE/SICE SII | C: gap G4 の主役。繊維入り材の試験片強度が歯で効かない実例 |
| 5 | Hriberšek ら「Technological optimization and fatigue evaluation of carbon reinforced polyamide 3D printed gears」 | 2024 | Heliyon | B: 「繊維入り」と「歯車疲労」を結ぶ。ただし懸念表明付きのため推奨値の導出には使わない |

深読みは `/paper-summary` を使わず、本節に直接書いた (Methodology 参照)。数値はすべて原典本文から取った。

### 5.1 Dennig2020 (FFF 歯車の VDI 2736 試験)

- 歯車: 中心距離 60mm、mn = 2mm、歯数 31 (印刷) と 29 (鋼)、歯幅 8mm (印刷) と 13mm (鋼)、圧力角 20°、平歯、転位なし (p5 表 1)。図 13 の注記では b = 7mm と書かれている (p11)。
- 印刷: すべて水平置き (歯車を寝かせる) でサポートなし。歯元に層間の引張がかからないようにした (p5 §3.2)。PLA は Craftbot+ で押出幅 0.4mm、層 0.2mm、外周 2、上下 4 層、充填 75%。ABS は Stratasys F270 で 0.3556mm、0.25mm、外周 2、上下 2 層、充填 22%。PA と PA+CF は Markforged Mark Two で 0.4mm、0.125mm、外周 2、上下 4 層、充填 37% (p6 表 3)。
- 試験: 鋼歯車と組ませた乾式運転、1500 rpm、各荷重 3 個以上、10⁶ または 2×10⁶ 回で打ち切り。印刷後 3 日以内に試験し、調湿していない (p11)。
- 結果 (p11 図 13、p12〜15):
  - PLA: 歯元温度が熱たわみ温度 63 °C を超えて歯が溶けた。σF ≈ 5〜11 MPa で ≈ 10⁴〜3×10⁴ 回 (図 13 の読み取り)。充填構造は寿命に影響しない (p12)。
  - ABS: すべての荷重で歯元折損。歯元温度は約 35 °C (p13)。図 13 の読み取りで、σF ≈ 13〜20 MPa で 1〜2×10⁴ 回、≈ 10〜17 MPa で 10⁵ 回前後、≈ 6〜9 MPa で 3〜5×10⁵ 回。外周を 2 から 4 に増やすと、σF ≈ 16 MPa で寿命が「わずかに」伸びた (p13)。
  - PA66: σF = 12 MPa で 2×10⁶ 回後も損傷なし。歯元温度は約 41 °C (p14 図 16)。
  - PA66+連続 CF: 2×10⁶ 回まで損傷なし。σF = 23 MPa を超えると約 60 °C に達して塑性変形し、≈ 31 MPa で破損した (p15)。
- 著者の結論: PLA は発熱で歯車に不適。ABS は低負荷なら使える。信頼できる値を得るには、歯元温度を一定にした S–N 曲線が要る (p15〜16)。
- 本調査への含意: 図 13 の値は図の目盛りからの読み取りで、±2 MPa 程度の誤差がある。PETG はこの試験に含まれない。

### 5.2 Martins2024 (PLA・PETG の S–N)

- 試験片: Anycubic Kobra、充填 100%、層 0.25mm、外殻 2mm、印刷速度 50 mm/s、ノズル 0.4mm (§2.1)。疲労試験片は 45° 方向に印刷 (要旨)。
- 引張強さ: PLA 55 MPa、PETG 37 MPa (平均、要旨と表 1)。
- 疲労 (R = 0.2、7 Hz、表 2・3):
  - PETG: σmax 20.4 MPa で 12,132〜15,483 回、15.3 MPa で 17,040〜41,250 回、10.2 MPa で 79,271〜193,323 回。
  - PLA: 32.4 MPa で 972〜2,040 回、24.3 MPa で 5,380〜6,880 回、16.2 MPa で 37,970〜48,893 回、10.8 MPa で 525,170 回 (1 本)。
  - PLA の Basquin 指数は −0.192 (§3.2 の考察)。
- 本調査での当てはめ (表 2・3 の全点に log σ = a + b log N を最小二乗で当てた。著者の当てはめではない):
  - PETG: 指数 −0.275。10⁴ 回 ≈ 21 MPa、10⁵ 回 ≈ 11 MPa、10⁶ 回 ≈ 5.9 MPa、10⁷ 回 ≈ 3.1 MPa。
  - PLA: 指数 −0.189。10⁴ 回 ≈ 22 MPa、10⁵ 回 ≈ 14 MPa、10⁶ 回 ≈ 9.2 MPa。
- 本調査への含意: 10⁶ 回以上は外挿である。荷重比 R = 0.2 の一方向引張で、歯元の繰返し (正逆転では引張と圧縮が入れ替わる) とは条件が違う。

### 5.3 Pogacnik2015 (規格値と実測値の比)

- 射出成形 POM (駆動) と PA6 (被動) の組。乾式、23 °C、600〜2305 rpm、0.30〜0.82 N·m (§3.4 表 2)。
- VDI 2736 で計算した歯元応力は 0.30 N·m で 17.6 MPa、0.42 で 24.6、0.59 で 34.6、0.82 で 48.1 MPa (p10 表 6)。
- 90% 生存率の寿命は 0.30 N·m・1646 rpm で 15.56×10⁶ 回、0.59 N·m・1176 rpm で 1.59×10⁶ 回 (p10 表 7)。
- 実測の許容値は VDI 2736 の PA66 値の 38〜72% (p10 表 7)。差は低荷重ほど大きい (0.30 N·m で 2.6 倍) (p10)。
- 本調査への含意: 規格の材料値を使う場合でも、実測で 0.4〜0.7 倍に下げる必要があった。印刷材ではこの比の情報すらない。

### 5.4 Yoshida2023 (印刷遊星 QDD)

- 減速機: 1 段遊星、減速比 8、歯数 太陽 9・遊星 27・内歯 63、平歯、m = 1.5、圧力角 20° (p2 表 I と本文)。理論定格出力トルク 37.6 N·m、最大 80 N·m (表 I)。
- 材料: Onyx (Markforged、短 CF 入りナイロン) と POTICON (チタン酸カリウム繊維入り)。曲げ強さは Onyx 71 MPa、POTICON 167 MPa、ABS 61 MPa、PLA 92 MPa (表 II)。
- 印刷: ノズル 0.4mm、押出温度 270 °C、層 0.1mm、外周 2、上下 4 層、充填 37% 三角 (表 IV)。
- 結果: 静的試験では 20 N·m 以下で効率 80% 超、測定範囲内で破損なし (p3)。正弦波運転 (最大 19.6 N·m) で POTICON の太陽歯車は約 2 時間で割れた。POM 太陽歯車の組では内歯車が約 46 時間で割れた (p5)。太陽歯車はフランジと歯の境界 (層の境界) で割れることが多く、表 II の強度に反して POTICON が Onyx より先に割れた。著者は層間接着の差で説明している (p4–5)。
- 本調査への含意: 繊維入り材の試験片強度は歯の強度の指標にならない。歯幅は論文に記載がなく、歯元応力は換算できなかった。

### 5.5 Hriberšek2024 (印刷 CF-PA 歯車、懸念表明付き)

- 歯車: z = 20、b = 6mm、m = 1、圧力角 20° (表 2)。Onyx を 0.4mm ノズル、充填 100% 三角で印刷 (表 3、§2)。相手は鋼歯車、1400 rpm (§2)。
- 寿命: 1.1 N·m で平均 6.27×10⁶ 回 (標準偏差 1.36×10⁶)、1.3 N·m で 2.05×10⁶ 回、1.5 N·m で 0.06×10⁶ 回 (§3)。1.5 N·m は過熱で樹脂が溶けたもので、疲労破壊ではないと著者は書いている (§3)。最適化前の歯車は歯元に充填不足があり、寿命が約 1/10 だった (§3)。
- 本調査での換算: 1.1 N·m のとき Ft = 110 N、本プロジェクトの Lewis 式で σ ≈ 54 MPa。これは同論文が測った Onyx の曲げ強さ 37.71 MPa (表 4) を上回り、Dennig2020 の PA66+CF (≈ 23 MPa で熱変形) とも合わない。
- 2025 年に Heliyon が懸念表明 (Expression of concern、doi:10.1016/j.heliyon.2025.e44032) を出している。理由は確認できなかった。このため推奨値の導出には使わない。

## 6. Paper Catalogue (RQ 別の要点)

### 6.1 歯車の試験 (RQ1)

- [[papers/Dennig2020_FffGearsFeasibility.pdf|Dennig 2020]] — FFF 歯車 (PLA・ABS・PA66・PA66+CF) の VDI 2736 疲労試験。§5.1。
- Pogacnik 2015 — 射出 PA6 歯車の加速試験と VDI 比較。§5.3。
- Hriberšek 2024 — 印刷 Onyx 歯車の寿命。懸念表明付き。§5.5。
- Tunalioglu & Agca 2022 — PLA・ABS・PETG の印刷歯車 (m = 6、歯数 17、歯幅 10mm、0.4mm ノズル、層 0.2mm、外周 4、充填 100%、印刷角 0°) を鋼歯車と組ませ、FZG 試験機で 1.5 N·m・900 rpm で運転した (表 2・3、§2.2)。摩耗は PETG が最少で、ABS より 30%、PLA より 40% 少ない (§3.1)。寿命は PETG が PLA より 22.3%、ABS より 37% 長い (§3.2)。寿命の絶対回数は図 7〜9 のみで、歯元応力は計算されていない。
- Nguyen 2025 — 中心距離 50mm、歯幅 12mm、圧力角 20° の印刷歯車の静的破壊試験 (p4 表 2)。印刷は層 0.2mm、55 mm/s、PETG は 240 °C (表 3)。m = 1.5 の平歯車 (33 歯) の破壊荷重は PLA 119 N、ASA 122 N、ABS 92 N、PETG 89 N (§3.3.3)。荷重はレバー端で測り、トルクは荷重 × 0.2 で表記されている (§2.2.3、§3.3)。ただし本文 PDF から単位とレバー長の単位 (mm か cm か) が抜けており、歯元応力に換算できない (未確認)。PETG の破壊荷重は PLA の約 0.75 倍で、ABS と同等だった。やまば歯 m = 1〜2.5 (PLA、充填 30%) の破壊荷重は 101〜151 N (表 6)。充填は 35% 以上、できれば 100% を推奨している (§3.3.2)。
- Kalani ら 2022 — 3D 印刷 PA12 平歯車 (m = 2、歯幅 10mm) の熱と摩耗。1400 rpm・2.4 N·m で 14,100 回で温度上昇により破損 (p6・p9)。印刷方式は FDM か SLS か本文で確認できていない (未確認)。
- Zhang, Mao, Leigh ら 2020 — 印刷条件の機械学習最適化で、ナイロン歯車の摩耗性能が 3 倍に改善 (要旨。本文未読)。
- Ciobanu ら 2025 — FDM PLA 歯車の摩耗。歯先が最も摩耗し、歯の内部は完全充填が必要と結論 (要旨。本文は数値未抽出)。
- Pisula ら 2021 — FDM・FFF 歯車列の摩耗を歯面形状で評価 (書誌のみ確認)。
- Zisopol ら 2022 — FDM 条件が PLA 平歯車の剛性に与える影響 (ローカル PDF。強度値は扱わない)。
- Dennig, Monn, Vodermayer 2019 — PEEK + 連続 CF の削り出し・射出歯車。VDI 2736 の材料データは PA66 と POM のみ、と記す (p2〜3)。印刷歯車ではないので比較参考。
- Pandian, Gautam, Senthilvelan 2022 — SLS PA12 歯車の層の向きと曲げ疲労 (検索結果の要約のみ。本文未確認、未検証)。

### 6.2 対象材の試験片 (RQ3・RQ4)

- Martins ら 2024 — PETG・PLA の S–N。§5.2。
- Dolzyk & Jung 2019 (J Fail Anal Prev 19(2):511–518) — PETG の 4 方向ラスタで UTS の 60〜90% の疲労試験 (二次情報のみ。未検証)。
- Kargar & Ghasemi-Ghalebahman 2023 — CF 30% 入り PLA の疲労。充填 50% では 45° ラスタ、75% では 0° ラスタが最長寿命。応力水準は UTS の 0.29〜0.52 倍 (要旨、§2)。
- CF-PETG の S–N (Int J Adv Manuf Technol, doi:10.1007/s00170-023-12107-4) — 4 応力水準の S–N を構築 (検索結果の要約のみ。著者と数値は未確認)。
- Valvez, Silva, Reis 2022 — 0.6mm 硬化鋼ノズル、ラスタ 0° の 3 点曲げ。最大曲げ強さは PETG 67.07 MPa (265 °C、30 mm/s、層 0.20mm、充填 100%)、PETG+CF 72.75 MPa (表 4・5)。全材料で充填率の影響が最大 (§3)。
- Economides ら 2024 — CF-PETG。充填 100% 縦ラスタの UTS 67.8 MPa が焼きなまし後 59.4 MPa (−12.3%) (§3)。分散寄与は充填パターン 22.01%、充填率 29.38%、焼きなまし 1.08% (§3)。
- Rihar & Hozdić 2026 — 線状充填 100% で引張強さ PETG 31.35 MPa、PETG+CF 38.90 MPa (要旨)。鉱物油に 7 日浸すと PETG は約 17% 低下、PETG+CF はほぼ不変 (要旨)。
- Hozdić & Hozdić 2025 — PLA+CF は充填 100% で引張強さ 33.4 MPa (要旨)。
- Spina 2026 — PET・PETG を 100 °C で焼きなますと引張強さ +8.1〜8.5%。190 °C では −22.7〜38.8% (§結果、表)。
- Zisopol ら 2026 — 焼きなまし PETG の最大曲げ強さ 84.78 MPa (層 0.20mm) (§3 表 2)。
- PETG・PETG-CF・PETG-KF の焼きなまし (Bhandari 等、Valvez2022 の文献紹介 §1): 炭素繊維入り PLA・PETG の引張強さが焼きなましで 2〜3 倍、と紹介されている (二次情報。Economides2024 とは逆向き)。

### 6.3 印刷減速機・アクチュエータ (RQ4)

- Yoshida ら 2023 — §5.4。
- Chi ら 2025 (Berkeley Humanoid Lite) — PLA 印刷のサイクロイド減速機。60 時間の耐久試験 (0.5 kg・0.5 m の振り子を −45°〜+90°、0.5 Hz) でバックラッシがわずかに増えた (§C Durability)。破損トルクの報告はない。
- Jenks ら 2025 (OpenSEA) — 印刷遊星 SEA。PLA の曲げ強さを 175 MPa と書いている (要旨の考察)。この値は一般的な PLA の値 (Yoshida2023 表 II で 92 MPa) より高く、根拠は確認できなかった。
- Bezzini ら 2026 — 印刷複合遊星 1:30 と印刷軸受。逆駆動トルク 0.42 ± 0.06 N·m (要旨)。歯の強度試験はない。
- Singh ら 2025 — 遊星の歯幅を Lewis 式で決める手法 (p3 式 9)。本プロジェクトの式の出典。

## 7. 設計への推奨

### 7.1 前提

- 応力の基準: プロジェクトのコードと同じ修正 Lewis 式 σ = Ft/(b·π·y·m)、y = 0.154 − 0.912/z、Kv = 1、ピーク入力 0.4 N·m、3 遊星で均等分担 (README.md 191 行)。
- 想定する使い方: 低速、ピークは間欠的。歯 1 枚が受けるピーク荷重の回数を **10⁵ 回** と仮定した (遠隔操作・軽量アーム関節で数年の使用を想定。この回数は本調査の仮定であり、文献値ではない)。
- 印刷: 歯車を寝かせて印刷、ノズル 0.4mm、歯の部分は実質充填 100% (外周 4 以上)。この条件を外れると下の値は使えない。

### 7.2 Lewis 値と VDI 値の換算 (本調査の概算)

文献の歯車データは VDI 2736 の σF で書かれている。Lewis 式との比は、本調査で次のように概算した。YFa・YSa は標準の歯形 (転位なし) の図表からの概略値で、規格原本とは照合していない (未確認)。

| 歯数 z | Lewis の係数 1/(π·y) | VDI の YFa·YSa·Yε (概算) | Lewis / VDI |
|---|---|---|---|
| 12 | 4.08 | ≈ 3.4 | ≈ 1.2 |
| 31 | 2.56 | ≈ 3.0 | ≈ 0.85 |

z = 12〜15 の太陽・遊星が最も応力が高く、そこでは Lewis が VDI より約 2 割大きく出る。つまり Lewis 基準の許容値を VDI 基準の試験値と同じ数字にすれば、小歯数の歯では約 2 割の余裕が上乗せされる。以下では換算せずに同じ数字として扱い、この余裕は安全率に含めない。

### 7.3 PETG

| 根拠 | 10⁵ 回で壊れる応力 | 出典 |
|---|---|---|
| PETG 試験片 (45° ラスタ、充填 100%、R = 0.2) | ≈ 11 MPa (当てはめ値)。実測点は 10.2 MPa で 0.8〜1.9×10⁵ 回 | Martins2024 表 2 |
| ABS 印刷歯車 (外周 2、充填 22%、σF) | ≈ 10〜17 MPa (図の読み取り) | Dennig2020 p11 図 13 |
| PETG と ABS の比 (歯車静的破壊) | PETG 89 N ≈ ABS 92 N | Nguyen2025 §3.3.3 |
| PETG と ABS の比 (歯車寿命) | PETG が ABS より 37% 長寿命 | Tunalioglu2022 §3.2 |

- 代表値 (10⁵ 回で壊れる応力): 上の下限側を取って **9 MPa**。
- 安全率: **1.5**。理由: (a) 低速で温度上昇による故障が支配的でない、(b) 小歯数では Lewis が約 2 割保守側、(c) 一方で、試験片から歯への換算が未確立で、層間割れ (Yoshida2023) と製造ばらつき (Hriberšek2024 で充填不足の歯車は寿命が約 1/10) がある。
- **推奨許容値: 9 / 1.5 = 6 MPa** (ピーク荷重 10⁵ 回)。
- 連続荷重 (10⁶ 回以上) に使う場合: 試験片の 10⁶ 回値 ≈ 5.9 MPa (Martins2024 表 2 からの当てはめ) を 1.5 で割り、**4 MPa**。
- 確度: **低〜中**。PETG 歯車そのものの S–N はなく、試験片と ABS 歯車からの推定である。

### 7.4 PETG-CF / PETG-GF

- 試験片の強度向上: 曲げ +8.5% (67.07 → 72.75 MPa、Valvez2022 表 4・5)、引張 +24% (31.35 → 38.90 MPa、Rihar2026 要旨)。
- 歯で強度が上がらない、または下がる根拠: 繊維入り材のカタログ強度が高くても先に割れた (Yoshida2023 p4–5)。焼きなましで CF-PETG の UTS が 12.3% 低下 (Economides2024 §3)。疲労試験片の散らばりは繊維と空隙で大きくなる (Kargar2023 の CF-PLA)。
- 歯車の S–N はない (gap G1)。GF 入り PETG の強度データは査読論文で見つからなかった。
- **推奨許容値: 6 MPa** (PETG と同じ。繊維による上乗せは認めない)。試験 (§7.7) で PETG に対して 1.2 倍以上の寿命が確認できたら、7 MPa まで上げてよい。
- 確度: **低**。

### 7.5 PLA-CF / PLA-GF

- 素の PLA 試験片は PETG より疲労に強い: 10⁵ 回 ≈ 14 MPa、10⁶ 回 ≈ 9.2 MPa (Martins2024 表 3 からの当てはめ)。
- しかし PLA の熱たわみ温度は 53〜63 °C (Tunalioglu2022 表 1、Dennig2020 p6 表 2) で、1500 rpm の試験では歯が溶けた (Dennig2020 p12)。ガラス転移温度は約 60 °C (Tunalioglu2022 表 1)。
- CF 入り PLA の引張強さは充填 100% で 33.4 MPa (Hozdić2025 要旨) で、素の PLA の 55 MPa (Martins2024 表 1) を下回る例がある。CF の追加で疲労寿命が 2〜21% 伸びた、という報告もある (Kargar2023 の考察)。
- 計算: PLA の 10⁵ 回値 14 MPa に、PETG で使った「試験片 → 歯車」の比 (9/11 ≈ 0.8) を掛けて約 11.5 MPa。安全率 1.5 で約 7.5 MPa。さらに温度とクリープ、繊維による脆化の不確かさで 0.8 倍し、**推奨許容値 6 MPa**。
- 条件: 歯車本体の温度を 40 °C 以下に保つこと (モータ発熱が伝わる配置では使わない)。この温度条件は、ガラス転移温度 60 °C に余裕を取った本調査の判断であり、文献値ではない。
- 確度: **低**。

### 7.6 現在の 10 MPa との比較

- 10 MPa は、PETG 試験片が 10⁵ 回前後で壊れる応力 (10.2 MPa で 0.8〜1.9×10⁵ 回、Martins2024 表 2) とほぼ同じで、安全率はほぼ 1 である。
- 許容値を 10 → 6 MPa に下げると、Lewis 式では歯幅 b は許容応力に反比例するので、歯幅は約 1.67 倍になる。
- 歯幅を増やせない場合の代わりの手段: やまば歯の重なりによる強度の割増し (現在のコードは割増しなし。README.md 191 行)、小歯車の正転位 (歯元が厚くなる)、歯元すみ肉半径を大きくする (ABS/CF-ABS 歯車ですみ肉 0.75mm が最良だった、という報告がある。Karupaiah ら 2024、BioResources 19(2):2796。本文未確認)。

### 7.7 検証の試験計画 (短縮版)

1. **静的破壊試験**: 実機と同じ太陽–遊星の組 (m 0.9 と 1.1、z 12〜15) を、PETG・PETG-CF・PLA-CF で各 5 個印刷する。遊星を固定し、太陽軸にトルクをかけて歯が折れるトルクを測る。Lewis 式で破壊応力に換算する。合格の目安: 破壊応力 ≥ 許容値 × 3 (= 18 MPa)。
2. **段階荷重試験** (Pogacnik2015 §3.3 の方法): 実機速度 (例: 入力 300 rpm 以下) で、Lewis 応力 4 MPa から 2×10⁴ 回ごとに 1 MPa ずつ上げ、破損した段を記録する。歯車本体の温度を赤外線で測る。
3. **寿命試験**: 段階試験の破損応力の 0.8・0.65・0.5 倍の 3 水準、各 3 個以上で、破損または 10⁶ 回まで回す。3 点から S–N を当てはめ、10⁵ 回値を許容値の代表値に置き換える (Dennig2020 p15〜16 が推奨する「歯元温度一定」の条件をできるだけ守る)。
4. **印刷条件の比較**: 外周 4 と 6、焼きなまし (PETG は 100 °C 前後。Spina2026) の有無を、手順 2 で比べる。
5. **判定**: 10⁵ 回値を 1.5 で割った値が §7.3〜7.5 の推奨値を上回れば、推奨値を測定値で置き換える。下回れば測定値に合わせて許容値を下げる。

## 8. Survey Methodology

### Search

- 情報源: WebSearch (一般)、Europe PMC REST API (オープンアクセス全文の取得と構造化検索)、ZHAW 機関リポジトリ (Dennig2020 の公開 PDF)、プロジェクト内の既存 PDF (`literature/papers/`)。
- 主な検索語: `3D printed PETG gear tooth bending fatigue VDI 2736`、`carbon fiber reinforced PETG FDM gear durability`、`FDM printed spur gear fatigue PLA PETG tooth root`、`TITLE:(gear) AND (3D print* OR FDM OR FFF)` (Europe PMC で 107 件)、`TITLE:(PETG) AND (fatigue OR fiber OR annealing)` (41 件)、`fatigue AND (carbon/glass fiber) AND (PLA OR PETG)` (18 件)。
- 雪だるま式の追跡: Tunalioglu2022 と Hriberšek2024 の引用から Dennig2020 を特定し、ZHAW リポジトリで本文を取得した。
- 取得できなかったもの: core.ac.uk (403)、TUM mediatum (ボット確認画面)、MDPI と PMC の直接取得 (403・reCAPTCHA)、Springer (認証へのリダイレクト)。いずれも回避はせず、Europe PMC の公開全文で代替できたものだけを使った。
- 採録数 27 本は focused の範囲 (20〜40) 内。歯車 × 疲労 × 対象材の組では直近の検索で新しい論文が出ず、飽和したと判断した。

### Map / Hub の確認

- 自動実行のため、Map Checkpoint と Hub 選定の承認は取っていない (automated run: presented but not confirmed)。
- Hub 深読みは `/paper-summary` を起動せず、本文から直接抽出して §5 に書いた (inline mode)。深読みノートを `literature/papers/{citekey}/` に分けて作る手順は省いた。
- 別エージェントによる深読みのクロス検証は行っていない。代わりに、§5 の数値はすべて原典本文 (PDF のテキストまたは Europe PMC 全文 XML) の該当箇所を grep で再確認した。

### Verify

- 書誌の存在確認: DOI のある論文は Europe PMC または出版社ページの記録で存在を確認した。`reference-verify` スキルは実行していない。
- 未検証の項目: Dolzyk & Jung 2019、Pandian ら 2022、CF-PETG の S–N (doi:10.1007/s00170-023-12107-4) の著者と数値、Karupaiah ら 2024 の数値、Pogacnik2015 の掲載先。これらは推奨値の導出に使っていない。
- 図からの読み取り値 (Dennig2020 図 13) は ±2 MPa 程度の誤差を含む。
- README と log の更新 (`literature/surveys/README.md`、`log.md`) は行っていない。

### 取得した PDF

- `literature/papers/Dennig2020_FffGearsFeasibility.pdf` (ZHAW 機関リポジトリの公開版。既存の運用に合わせてローカルに置いた。再配布はしない)。

## 9. 参考文献

1. Dennig, H.-J., Zumofen, L., Kirchheim, A. (2020). Feasibility investigation of gears manufactured by fused filament fabrication. In: 2nd Int. Conf. on Additive Manufacturing for Products and Applications (AMPA 2020), Springer, pp. 304–320. ISBN 978-3-030-54333-4. https://doi.org/10.21256/zhaw-20460
2. Pogačnik, A., Tavčar, J. (2015). Accelerated testing method for polymer gears — Review of gear body temperature calculation according to VDI 2736. 会議論文 (掲載先未確認)。ローカル: `literature/papers/Pogacnik2015_PolymerGearAcceleratedTesting.pdf`
3. Hriberšek, M., Kulovec, S., Ikram, A., Kern, M., Kastelic, L., Pušavec, F. (2024). Technological optimization and fatigue evaluation of carbon reinforced polyamide 3D printed gears. Heliyon 10, e34037. https://doi.org/10.1016/j.heliyon.2024.e34037 (Expression of concern: https://doi.org/10.1016/j.heliyon.2025.e44032)
4. Tunalioglu, M.S., Agca, B.V. (2022). Wear and Service Life of 3-D Printed Polymeric Gears. Polymers 14(10), 2064. https://doi.org/10.3390/polym14102064
5. Nguyen Von Dim (2025). Experimental investigation of 3D printed gears. TNU Journal of Science and Technology 230(10), 463–470. ローカル: `literature/papers/Nguyen2025_3dPrintedGearsExperimental.pdf`
6. Yoshida, T., Endo, G., Okubo, A., Nabae, H. (2023). Experimental Evaluation of a Quasi-direct-drive Actuator with a 3D-printed Planetary Gear Reducer. 2023 IEEE/SICE SII. https://doi.org/10.1109/SII55687.2023.10039274
7. Martins, R.F., Branco, R., Martins, M., Macek, W., Marciniak, Z., Silva, R., Trindade, D., Moura, C., Franco, M., Malça, C. (2024). Mechanical Properties of Additively Manufactured Polymeric Materials—PLA and PETG—For Biomechanical Applications. Polymers 16(13), 1868. https://doi.org/10.3390/polym16131868
8. Valvez, S., Silva, A.P., Reis, P.N.B. (2022). Optimization of Printing Parameters to Maximize the Mechanical Properties of 3D-Printed PETG-Based Parts. Polymers 14(13), 2564. https://doi.org/10.3390/polym14132564
9. Economides, A.L., Islam, M.N., Baxevanakis, K.P. ほか (2024). Additively Manufactured Carbon Fibre PETG Composites: Effect of Print Parameters on Mechanical Properties. Polymers 16, 3336. https://doi.org/10.3390/polym16233336
10. Rihar, L., Hozdić, E. (2026). Mechanical and Microstructural Response of FDM-Printed PETG and PETG+CF to Variable Infill Architecture and Lubricant Exposure. Polymers 18(5), 654. https://doi.org/10.3390/polym18050654
11. Hozdić, E., Hozdić, E. (2025). Mechanical Parameters and Microstructural Evolution of FDM-Printed PLA and PLA+CF Under Variable Infill Architecture and Lubricant Exposure. Polymers 18(1), 72. https://doi.org/10.3390/polym18010072
12. Spina, R. (2026). Annealing of PET/PETG specimens made with material extrusion. Scientific Reports. https://doi.org/10.1038/s41598-026-50018-y
13. Zisopol, D.G., Minescu, M., Iacob, D. ほか (2026). Influence of Post-Processing by Annealing on Additive Manufacturing Productivity and Mechanical Properties in Three-Point Bending of Additively Manufactured PETG and rPETG. Polymers. https://doi.org/10.3390/polym18192348
14. Kargar, E., Ghasemi-Ghalebahman, A. (2023). Experimental investigation on fatigue life and tensile strength of carbon fiber-reinforced PLA composites based on fused deposition modeling. Scientific Reports 13. https://doi.org/10.1038/s41598-023-45046-x
15. Kalani, A., Vadher, J., Sharma, S., Jani, R. (2022). Investigation of Thermal and Wear Behaviour of 3D Printed PA-12 Nylon Polymer Spur Gears. El-Cezerî J. Sci. Eng. 9(3), 1121–1135. https://doi.org/10.31202/ecjse.1085328
16. Zhang, Y., Mao, K., Leigh, S., Shah, A., Chao, Z., Ma, G. (2020). A parametric study of 3D printed polymer gears. Int. J. Adv. Manuf. Technol. 107. https://doi.org/10.1007/s00170-020-05270-5
17. Ciobanu, R., Arhip, G., Donțu, O., Rizescu, C.I. ほか (2025). Additive Manufacturing Meets Gear Mechanics: Understanding Abrasive Wear Evolution in FDM-Printed Gears. Polymers 17(13), 1810. https://doi.org/10.3390/polym17131810
18. Pisula, J., Budzik, G., Turek, P. ほか (2021). An Analysis of Polymer Gear Wear in a Spur Gear Train Made Using FDM and FFF Methods Based on Tooth Surface Topography Assessment. Polymers 13(10), 1649. https://doi.org/10.3390/polym13101649
19. Zisopol, D.G., Iacob, D.V., Portoaca, A.-I. (2022). A Theoretical-Experimental Study of the Influence of FDM Parameters on PLA Spur Gear Stiffness. Eng. Technol. Appl. Sci. Res. 12(5), 9329–9335. ローカル: `literature/papers/Zisopol2022_PlaSpurGearFdmParameters.pdf`
20. Dennig, H.-J., Monn, S., Vodermayer, A. (2019). Thermoplastic high performance composite gears. 3rd Int. Conf. on High Performance Plastic Gears. https://digitalcollection.zhaw.ch/bitstream/11475/18362/1/3rd%20International%20Conference%20on%20Gears%20-%20DENNIG%20-%20final.pdf
21. Pandian, A.K., Gautam, S.S., Senthilvelan, S. (2022). SLS Nylon 12 歯車の層の向きと曲げ疲労 (題名・掲載先 未検証)。
22. Dolzyk, G., Jung, S. (2019). Tensile and Fatigue Analysis of 3D-Printed Polyethylene Terephthalate Glycol. J. Fail. Anal. Prev. 19(2), 511–518 (未検証。二次情報のみ)。
23. CF-PETG の S–N 曲線 (2023). Int. J. Adv. Manuf. Technol. https://doi.org/10.1007/s00170-023-12107-4 (著者・題名 未検証)。
24. Chi, Y., Liao, Q., Long, J., Huang, X., Shao, S., Nikolic, B., Li, Z., Sreenath, K. (2025). Demonstrating Berkeley Humanoid Lite. arXiv:2504.17249.
25. Jenks, B., Levan, H., Stefanovic, F. ほか (2025). OpenSEA: a 3D printed planetary gear series elastic actuator for a compliant elbow joint exoskeleton. Front. Robot. AI. https://doi.org/10.3389/frobt.2025.1528266
26. Bezzini, R., Bassani, G., Avizzano, C.A. ほか (2026). Design and Examination of Compound Planetary Gearboxes to Investigate the Impact of 3D-Printed Embedded Bearings for Assistive Exoskeletons Actuators. J. Mech. Des. 148(11), 113501. https://doi.org/10.1115/1.4071656
27. Singh, A., Kapa, D., Chedda, P., Kolathaya, S.N.Y. (2025). Comparison between External and Internal Single Stage Planetary gearbox actuators for legged robots. arXiv:2506.16356.
28. Karupaiah ほか (2024). Performance Evaluation of 3D-Printed ABS and Carbon Fiber-reinforced ABS Polymeric Spur Gears. BioResources 19(2), 2796 (本文未確認)。

## 略語一覧

| 略語 | 完全名 | 初出 |
|---|---|---|
| FDM / FFF | Fused Deposition Modeling / Fused Filament Fabrication (熱溶解積層) | §0 |
| CF / GF | Carbon Fiber / Glass Fiber (炭素繊維 / ガラス繊維) | §0 |
| PETG | Polyethylene Terephthalate Glycol-modified | §0 |
| PLA | Polylactic Acid (ポリ乳酸) | §0 |
| PA | Polyamide (ナイロン) | §1 |
| S–N | 応力–繰返し数曲線 (Wöhler 曲線) | §1 |
| UTS | Ultimate Tensile Strength (引張強さ) | §6.2 |
| QDD | Quasi-Direct-Drive (準ダイレクトドライブ) | §5 |
| SLS | Selective Laser Sintering (粉末焼結) | §0 |
| VDI 2736 | ドイツ技術者協会の樹脂歯車設計指針 | §0 |
