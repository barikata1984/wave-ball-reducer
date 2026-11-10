# GM5208-24 用 1 段遊星減速機 (FDM、やまば歯) の設計根拠 (2026-11-10)

目的: iFlight GM5208-24 を QDD 化する 1 段遊星減速機 (減速比 6〜8、PETG/PLA、0.4 mm ノズル) を build123d でパラメトリック設計するための式と数値を、出典付きでまとめる。前段の方式選定は `2026-11-10_printable-qdd-reducer.md` を参照 (重複は省く)。

凡例: [事実] は出典とページ付き。[推測] は筆者の推論・概算。ページ番号は PDF のページ (1 始まり)。

## 1. 文献一覧

| 分野 | 書誌 | URL/DOI | ローカル PDF |
|---|---|---|---|
| A | KHK (小原歯車工業), "Internal Ring Gears" catalog 2025 (SI/SIR 章、Application Hints を含む) | https://khkgears.us/media/404afki4/catalog-2025-internal-ring-gears.pdf | literature/papers/Khk2025_InternalRingGearsCatalog.pdf |
| A | KHK, Gear Technical Reference "Involute gear profile" (Web) | https://khkgears.net/new/gear_knowledge/gear_technical_reference/involute_gear_profile.html | 未入手: Web ページのみ (PDF なし) |
| A | KHK, "Addendum modification and undercut" | https://www.khkgears.us/media/1284/addendum-modification-and-undercut.pdf | 未入手: ダウンロードすると 0 バイト |
| A/B/C | Singh, Kapa, Chedda, Kolathaya, "Comparison between External and Internal Single Stage Planetary gearbox actuators for legged robots," arXiv:2506.16356, 2025 | https://arxiv.org/abs/2506.16356 | literature/papers/Singh2025_ExternalInternalPlanetary.pdf |
| B/D | Nguyen Von Dim, "Experimental investigation of 3D printed gears," TNU J. Sci. Tech. 230(10):463–470, 2025 (本文ベトナム語) | https://jst.tnu.edu.vn/jst/article/download/12714/pdf | literature/papers/Nguyen2025_3dPrintedGearsExperimental.pdf |
| B | Zisopol et al., "A Theoretical-Experimental Study of the Influence of FDM Parameters on PLA Spur Gears Stiffness," ETASR 12(5):9329–9335, 2022 | https://distantreader.org/stacks/journals/etasr/etasr-5183.pdf | literature/papers/Zisopol2022_PlaSpurGearFdmParameters.pdf |
| B | VDI 2736 (Thermoplastic gear wheels) | 規格 (有料) | 未入手: 有料規格。規格の一部の値は下の Pogacnik ら経由で確認 |
| B | Pogacnik, Tavcar, "Accelerated testing method for polymer gears — Review of gear body temperature calculation according to VDI 2736" (VDI 会議の原稿。PDF に誌名・年・DOI の記載なし。ファイル名と本文 "In the last 2 years, the new VDI 2736" から 2015 年の VDI 歯車会議と推定 [推測]、会議名は未検証) | DOI なし (未検証)。関連する査読論文は Pogacnik & Tavcar, Materials and Design 65 (2015) 961–973 | literature/papers/Pogacnik2015_PolymerGearAcceleratedTesting.pdf (元 Gear_VDI_2015_Manuscript_v6_errata_cor.pdf) |
| C | Chi, Liao, Long, Huang, Shao, Nikolic, Li, Sreenath, "Demonstrating Berkeley Humanoid Lite," RSS 2025, arXiv:2504.17249 | https://arxiv.org/abs/2504.17249 | literature/papers/Chi2025_BerkeleyHumanoidLite.pdf |
| C | Singh, Kapa, Joshi, Kolathaya, "COMPAct: Computational Optimization and Automated Modular design of Planetary Actuators," arXiv:2510.07197v2, 2026 (ICRA 2026) | https://arxiv.org/abs/2510.07197 | literature/papers/Singh2025_COMPAct.pdf |
| C | Katz, Di Carlo, Kim, "Mini Cheetah: A Platform for Pushing the Limits of Dynamic Quadruped Control," ICRA 2019, pp. 6295–6301 (PDF 1 枚目で題・著者・会議を確認。DOI は PDF になく IEEE Xplore 8793865 で確認) | DOI 10.1109/ICRA.2019.8793865 | literature/papers/Katz2019_MiniCheetah.pdf |
| C | Yoshida, Endo, Okubo, Nabae, "Experimental Evaluation of a Quasi-direct-drive Actuator with a 3D-printed Planetary Gear Reducer," IEEE/SICE SII 2023 (PDF 1 枚目で確認) | DOI 10.1109/SII55687.2023.10039274 | literature/papers/Yoshida2023_3dPrintedPlanetaryQdd.pdf |
| C | Bezzini, Bassani, Avizzano, **Filippeschi**, "Design and Examination of Compound Planetary Gearboxes to Investigate the Impact of 3D-Printed Embedded Bearings for Assistive Exoskeletons Actuators," J. Mech. Des. 148(11):113501, Nov. 2026 (online 2026-04-27)。想定どおりの論文。著者は 4 名 | DOI 10.1115/1.4071656 | literature/papers/Bezzini2026_CompoundPlanetaryPrintedBearings.pdf (元 md-25-1816.pdf) |
| D | Polymers 17(3):416, 2025, "Geometric Accuracy and Dimensional Precision in 3D Printing-Based Gear Manufacturing" | DOI 10.3390/polym17030416 (PMC11820634) | 未入手 (PDF): PMC が自動取得を拒否。Europe PMC の全文 XML で本文を確認 (ページ番号なし) |
| E | iFlight GM5208-24 製品ページ | https://shop.iflight.com/ipower-motor-gm5208-24-brushless-gimbal-motor-pro1347 | 対象外 (Web) |
| B | 材料 TDS (Polymaker PolyLite PETG/PLA、eSUN PETG) | polymaker.com / esun3d.com | 未入手: サイトが PDF でなく HTML を返す |

## 2. 設計に使う式と数値

記号: Zs 太陽、Zp 遊星、Zr 内歯車の歯数、N 遊星数、m モジュール、α 圧力角。

### 2.1 遊星の歯数条件
- [事実] 同心条件 Zr = Zs + 2Zp。組立条件 (Zs + Zr)/N = 整数。隣接条件 Zp + 2 < (Zs + Zp)·sin(180°/N)。KHK catalog p3 (Application Hints ③)。Singh2025_ExternalInternalPlanetary p3 も同じ同心条件と組立条件を示す。
- [事実] 隣接条件の別形 (遊星どうしの最小すき間 δp を明示): 2m(Zs+Zp)·sin(π/N) − 2mZp ≥ δp。同論文では δp = 5 mm。Singh2025_ExternalInternalPlanetary p4 式 (14)。
- [事実] 減速比 (内歯車固定、太陽入力、キャリア出力) は 1 + Zr/Zs。効率は η = (Zs + ηsp·ηpr·Zr)/(Zs + Zr)。Singh2025_COMPAct p2 式 (2)。
- [推測] 候補 (N = 3、α = 20°):
  - 6:1 は Zs 12、Zp 24、Zr 60。組立 72/3 = 24、隣接 26 < 36·0.866 = 31.2 で成立。
  - 8:1 は Zs 12、Zp 30、Zr 72 (比 7.0) または Zs 12、Zp 36、Zr 84 (比 8.0)。後者は組立 96/3 = 32、隣接 38 < 41.6 で成立。m = 0.8 なら内歯車ピッチ円は 67.2 mm。

### 2.2 最小歯数・転位・内歯車干渉
- [事実] α = 20° の標準歯車で、切り下げの起きない最小歯数は 17。正の転位で切り下げを防げる。KHK Gear Technical Reference (Web、ページ番号なし)。
- [推測] 理論式 z_min = 2/sin²α = 17.1。切り下げを防ぐ転位係数は x ≥ 1 − z·sin²α/2 で、z = 12 なら x ≥ 0.30、z = 10 なら x ≥ 0.42。式は標準的な教科書式で、KHK の該当 PDF は入手できず未照合。
- [事実] 歯数 60 の内歯車に噛み合う小歯車の許容範囲 (KHK 在庫品の条件): インボリュート干渉の下限 21、トロコイド干渉の上限 51、トリミング干渉の上限 43。歯数 80 では下限 20、上限 72 と 64。KHK catalog p2 の表。
  - 原因の整理: インボリュート干渉は小歯車の歯数が少なすぎるとき、トロコイド干渉とトリミング干渉は両歯車の歯数差が小さすぎるときに起きる。同 p2。
- [推測] 上の表は KHK 製品の歯形の前提での値。転位なしの標準歯形で作るなら、遊星 Zp ≥ 21 (Zr 60 のとき) を目安にし、最終確認は build123d で歯形どうしの干渉を幾何チェックする。
- [事実] 比較研究では、切り下げを避ける目的で太陽・遊星とも Zmin = 20、モジュール 0.5〜1.2 mm に制約している。Singh2025_ExternalInternalPlanetary p4。

### 2.3 プラスチック歯車の強度
- [事実] 修正 Lewis 式 σ = Ft/(b·y·Kv·P) (P は円ピッチ、y は Lewis 形状係数、Kv は速度係数)。歯幅 b はこの式から決める。Singh2025_ExternalInternalPlanetary p3 式 (9)。COMPAct も同じく修正 Lewis 式で歯幅を決める (Singh2025_COMPAct p4)。
- [事実] 太陽の許容伝達トルク T1 = Ts·Zp·η (Ts は太陽と遊星 1 組の許容トルク、Zp は遊星の数)。KHK catalog p5 式 (1)。
- [事実] 印刷歯車の破壊荷重 (PLA、層厚 0.2 mm、インフィル 30 %、やまば歯):

  | m | 破壊荷重 | トルク |
  |---|---|---|
  | 1.0 | 101 N | 20.2 |
  | 1.25 | 129 N | 25.8 |
  | 1.5 | 110 N | 22.0 |
  | 2.0 | 151 N | 30.2 |
  | 2.5 | 124 N | 24.8 |

  Nguyen2025 p7 表 6。トルクの単位は PDF で文字化けして読めない (荷重 × 0.2 と一致するので、腕 0.2 m の N·m と推定 [推測])。
- [事実] 材料比較 (平歯、m = 1.5): PLA 119 N、ASA 122 N、ABS 約 92 N、PETG 89 N で破壊。歯形比較 (m = 2.5): 平歯 192 N、はすば 210 N、やまば 149 N。インフィルは 35 % 以上を推奨。Nguyen2025 p7。
- [事実] VDI 2736 と JIS B 1759 の許容曲げ応力は、歯車試験機での実測値に寿命係数を掛けて求める。VDI のデータは POM、PA66、PET、PE の 4 材料だけで、PLA と PETG はない (Gear Solutions 記事 https://gearsolutions.com/?p=53751、二次資料)。
- [事実] VDI 2736 の許容歯元応力 σFE (PA66) は 46〜52 MPa (61〜78 °C、寿命 1.6〜15.6×10^6 回の条件)。同じ条件で射出成形 PA6 (相手 POM、m 1、z 20、歯幅 6 mm、無潤滑) を歯車試験した実測の許容値 (90 % 生存) は 17.6〜34.6 MPa で、VDI 値の 38〜72 %。差は低荷重・長寿命ほど大きい。Pogacnik p10 表 7。著者は「VDI 2736 の許容値は楽観的すぎる可能性」を挙げる (p9–10, p11)。
- [事実] 歯元応力の換算 (VDI 2736 式、m 1、z 20、b 6 mm): 0.30 N·m → 17.6 MPa、0.42 → 24.6、0.59 → 34.6、0.82 → 48.1 MPa。Pogacnik p9 表 6。歯元温度の式 (VDI 2736) は p8 式 (1) だが画像で、本文抽出できず未転記。POM/PA の標準係数は μ 0.18、kθ,Fuß 2100、c 0.75。実測に合わせると μ 0.36、kθ,Fuß 1505、c 0.36 (標準偏差 20 °C → 4 °C)。p8–9 表 5。
- [事実] 破損モードは歯元疲労と熱による溶融の 2 種。溶融は過負荷と温度上昇で起き、多くは試験開始 1 時間以内。それを越えて温度が安定すると疲労が主。摩耗は材料の組合せ (POM/POM) 次第。同材 PA6/PA6 は摩擦係数が高く 0.45 N·m で熱破損。Pogacnik p1, p6, p10。
- [事実] 歯車試験の段階荷重法: 開始時の歯元応力 15〜20 MPa (PA6/POM)、1 段あたり 20 % 増し、各段 2×10^5 回以上。Pogacnik p3。
- [事実] 熱膨張を見込み、中心距離を 20.00 mm (標準) から 20.05 mm に広げて運転。Pogacnik p5 表 1 の前文。
- [事実] 印刷遊星 (Onyx / POTICON、インフィル 37 % 三角、シェル 2、上下 4 層、0.4 mm ノズル、積層 0.1 mm、270 °C): 太陽歯車は「歯とフランジの境界 (積層面)」で破断。POTICON はカタログ強度が Onyx の 2 倍以上なのに太陽が弱く、19.6 N·m 正弦負荷で約 2 時間で破断。POM 切削太陽の組では印刷内歯車が約 46 時間で破断し、歯の内部はインフィルがほとんど入っていなかった。著者の結論は「耐久性はフィラメント物性よりスライサの経路と内部充填が支配的」。Yoshida p4–5、表 IV (p5)。
- [推測] 上の 2 文献から、FDM PETG の許容曲げ応力は「材料の引張強さ × 係数」では決まらない。射出成形品でも歯車試験の許容値は規格値の 4〜7 割、印刷品は歯の内部充填と積層面で更に落ちる。PLA/PETG の印刷歯車で許容曲げ応力を歯車試験で実測した一次資料は依然見つかっていない。概算例: 入力 0.25 N·m、Zs 12、m 1.0、N 3 なら、1 噛み合いあたり Ft ≈ 0.25/0.006/3 ≈ 14 N。b = 10 mm、Y ≈ 0.25 で σ ≈ 5.6 MPa。PETG の引張強さ約 47 MPa (第三者の販売ページの値、未検証) に対して余裕は大きい。積層間の強度はその半分程度とみて、許容 10 MPa 程度を仮の設計値にし、試作で破壊試験をする。

### 2.4 バックラッシ・クリアランス・寸法精度
- [事実] 印刷歯車対の両方向合計バックラッシ (ISO TR 10064-2 に準拠した測定): 平歯 2.6、はすば 1.6、やまば 1.9。いずれも「距離 10」の位置で測定しており、単位は文字化けで読めない。Nguyen2025 p6 表 4。
- [事実] FDM 歯車 (m = 2) の寸法偏差は ±0.1 mm 以内に収まり、小さい部品ほど熱収縮と異方性の影響を受ける。Polymers 17(3):416 (全文 XML、ページなし)。
- [事実] 内歯車は、中心距離 a = (d2 − d1)/2 で組んだときに適正なバックラッシが出るよう設計されている。KHK catalog p3。
- [事実] Yoshida の印刷遊星はバックラッシ量の設計値・実測値を示していない。剛性低下の理由に「造形精度によるバックラッシ」を挙げるのみ。Yoshida p4。
- [事実] Bezzini の印刷複合遊星 (PLA) のギアプレイ (出力の遊び角): 軸受版 0.114 ± 0.076°、印刷軸受版 0.195 ± 0.084°、先行 CO30 0.53 ± 0.21°。他の印刷減速機の最良値は 0.23〜0.7°。Bezzini p7 表 3。歯厚減らし・中心距離の設計値は本文に記載なし。印刷軸受のローラ径は 3 mm を 3.05 mm として軌道を設計 (片側 0.025 mm 相当 [推測])。Bezzini p3。
- [推測] 推奨クリアランスの一次資料はない。既存リポジトリのクリアランス (片側 0.05 mm、コミット 0fb9d00) を初期値とし、中心距離を +0.05〜0.15 mm ずらした試験片を印刷して決める。

### 2.5 FDM の最小モジュール・歯幅
- [事実] 比較研究の製造制約は m = 0.5〜1.2 mm。Singh2025_ExternalInternalPlanetary p4。
- [事実] PLA やまば歯の破壊荷重は m = 1.0 で 101 N。Nguyen2025 p7。
- [事実] COMPAct の印刷条件はインフィル 100 %、Bambu Lab X1C。Singh2025_COMPAct p6。
- 経験則 (m 0.65 の印刷実例、歯幅 m の 12〜20 倍) は前段の調査 §5 を参照。
- [推測] 0.4 mm ノズルでは m 0.8〜1.0 を下限とする。歯幅は 2.3 の概算から 8〜12 mm で足りる。

### 2.6 事例の実測値
| 事例 | 減速機 | 効率 | バックラッシ | 剛性 | 出典 |
|---|---|---|---|---|---|
| Berkeley Humanoid Lite 6512 | 印刷サイクロイド、PLA | 機械効率 約 90 % | 最大 0.0229 rad (SD 0.0042)、6 台 | 319.49 N·m/rad | Chi2025 p6 (効率と剛性)、p7 (バックラッシ) |
| COMPAct SSPG (7.2:1) | 印刷 1 段遊星 | 機械効率 60〜80 % (モデルの予測は 96 %) | 無負荷 0.010 rad (0.59°) | 242.7 N·m/rad | Singh2025_COMPAct p1, p6–7 |
| COMPAct CPG (14:1) | 印刷複合遊星 | 約 60 % | 0.046 rad (2.6°) | 201.6 N·m/rad | 同上 |
| Bezzini ら C3D | 印刷複合遊星 30:1、PLA (Ultimaker S7、Cura)、m 1.2、α 25°、遊星 4、歯数 Zsa 16 / Zpa 16 / Zra 48 / Zpb 13 / Zrb 45 | — (摩擦: 起動 27.59 mN·m) | 0.195 ± 0.084° | 653.68 ± 64.91 N·m/rad | Bezzini p3 表 1、p4、p7 表 3、p8 表 4、p9 表 5。逆駆動トルク 0.42 ± 0.063 N·m (p7 表 3)。軸受版 CFB は 1.18 ± 0.019 N·m、0.114°、971.52 N·m/rad |
| Yoshida ら (改良 OpenTorque) | 印刷 1 段遊星 8:1、平歯 (はすばから変更)、m 1.5、α 20°、Zs 9 / Zp 27 / Zr 63、Onyx・POTICON | 静トルク 20 N·m 以下で 80 % 超、それ以上で低下 | 値なし | 数値なし (図のみ) | Yoshida p1 表 I、p2、p3。定格 37.6 N·m・最大 80 N·m は理論値 |
| Mini Cheetah | 1 段遊星 6:1、Misumi の市販金属歯車、遊星はステータ内に配置、ニードル軸受 | 値なし | 値なし | 値なし | Katz p3。最大 17 N·m、連続 6.9 N·m (p4 表 I)。膝はさらにベルト 1.55:1 (p2)。伝達部の許容ピークトルクはモータ最大の約 4 倍 (p4)。歯数・モジュールは記載なし |

- [事実] Yoshida は軸方向力を避けるため、遊星をはすばから平歯に変えた。遊星は 2 個の玉軸受で支持して傾きを抑えた。Yoshida p2。
- [事実] Yoshida の 46 時間耐久試験では、内歯車表面は約 55 °C で熱変形温度未満。Yoshida p4。

- [事実] 60 時間の耐久試験で、効率はほぼ維持され、バックラッシは摩耗で徐々に増えた。Chi2025 p6–7。
- [事実] 逆駆動トルクを測った印刷 1 段遊星の一次資料は、今回入手した 4 件にもない。最も近いのは Bezzini の 30:1 複合遊星 (0.42 N·m、出力側)。Bezzini p7 は 1 段遊星 QDD の 0.97 N·m の例も挙げる (出典番号のみ確認、元論文は未確認)。

## 3. GM5208-24 の取付寸法

- 公式: **不明**。製品ページに取付穴の記載も図面もない (2026-11-10 確認)。
- 参考 (別機種、確度は低い): iPower GBM5208H は上面 2 穴 12 mm、3 穴 17 mm (中心間)、底面 3 穴 32 mm、4 穴 25 mm (中心間)。出典は販売店ページ (https://www.onedrone.com/store/ipower-gbm5208h-200t-hollow-shaft.html) と manuals.online の GBM5208-150 マニュアル。中空軸の寸法が GM5208-24 と異なるので、流用できない可能性がある。
- [推測] 実物をノギスで測るまで、build123d の取付穴はパラメータ (PCD、穴数、ねじ径) として未確定のまま残す。

## 4. 未入手リスト
1. KHK "Addendum modification and undercut" PDF: 0 バイトで返る
2. VDI 2736: 有料規格
3. (2026-11-10 入手済み: Katz 2019、Yoshida 2023、Bezzini 2026、Pogacnik の VDI 会議原稿)
4. VDI 2736 の歯元温度の式本体: Pogacnik p8 では画像で抽出できず。規格本体は有料
5. Polymers 17(3):416 の PDF: PMC が自動取得を拒否 (XML 全文で代用)
6. PETG/PLA 材料 TDS (Polymaker、eSUN、Prusament): サイトが HTML を返す
7. GM5208-24 の取付図: 公開資料が見つからない
