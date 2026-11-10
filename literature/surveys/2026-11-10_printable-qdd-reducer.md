# 小型 FDM プリンタで印刷できる QDD 用減速機の調査 (2026-11-10)

目的: FDM 3D プリンタで印刷でき、アウターロータ BLDC を QDD (quasi-direct drive: 低減速比で逆駆動性と力の透明性を保つアクチュエータ) にする減速機の方式を選ぶ。前提機材は iFlight GM5208-24 と Creality Hi。

凡例: [事実] は出典付き。[推測] は筆者の推論で、出典はない。[未検証] は一次資料で確認できなかった値。

## 1. 前提機材の仕様

### iFlight iPower GM5208-24 (出典: https://shop.iflight.com/ipower-motor-gm5208-24-brushless-gimbal-motor-pro1347)
| 項目 | 値 |
|---|---|
| 外径 | Φ59.5 ±0.05 mm |
| 高さ | 31.2 ±0.2 mm |
| 軸 | 中空軸。外径 Φ15、内径 Φ12.6 mm |
| 構成 | 24N/22P (11 極対) |
| 巻線抵抗 | 13.7 Ω ±5% |
| KV | 記載なし。無負荷回転数 396〜436 rpm @20 V から換算すると約 20〜22 rpm/V [推測: 換算値] |
| トルク | 負荷トルク 1800〜2500 g·cm (約 0.18〜0.25 N·m)。最大トルクの記載なし |
| 電圧・電力 | 3〜5S、最大 40 W、負荷電流 1 A @20 V |
| 質量 | 204 g |
| 取付穴 | 公式ページに記載なし。実物の実測か図面の入手が必要 |

[推測] 13.7 Ω の高抵抗巻線なので、24 V でも電流は約 1.75 A にとどまる。連続トルクは 0.2 N·m 前後、短時間でも 0.3〜0.4 N·m 程度が上限とみられる。

### Creality Hi (出典: https://www.creality.com/support/creality-hi-combo)
| 項目 | 値 |
|---|---|
| 造形範囲 | 260 × 260 × 300 mm (依頼時の想定 180 mm 角より大きい) |
| ノズル | 0.4 mm 標準。0.6 / 0.8 mm に対応。トライメタル製 |
| 最高温度 | ノズル 300 °C、ベッド 100 °C |
| 公式の対応材料 | PLA, TPU 95A, BVOH, PLA-CF, PVA, PETG, ABS |
| 筐体 | オープンフレーム (囲いなし) |

[推測] 筐体がないため、PA (ナイロン) や ABS では反りが出やすい。PA-CF は小売店の記載 (https://www.3djake.com/creality-3d-printers-spare-parts/hi-combo) にはあるが、公式の一覧にはない。

## 2. QDD の減速方式と逆駆動性の要因

- [事実] MIT Cheetah 系は、高トルク密度の大径モーターに低減速比を組み合わせ、逆駆動性を確保した。逆駆動性の指標として impact mitigation factor を提案している。Seok, Wensing, Wang, Otten, Lang, Kim, "Proprioceptive Actuator Design in the MIT Cheetah," IEEE T-RO 33(3):509–522, 2017. https://dspace.mit.edu/handle/1721.1/119863
- [事実] Mini Cheetah のアクチュエータは遊星歯車を内蔵し、ピークトルクは 17 N·m、連続トルクは 6.9 N·m。Katz, Di Carlo, Kim, "Mini Cheetah: A Platform for Pushing the Limits of Dynamic Quadruped Control," ICRA 2019, pp. 6295–6301. https://dspace.mit.edu/handle/1721.1/127165 ; 修論は Katz, "A Low Cost Modular Actuator for Dynamic Robots," MIT 2018. https://dspace.mit.edu/handle/1721.1/118671
- [未検証] 減速比 6:1 は広く引用される値だが、今回は一次資料の本文で確認できなかった (DSpace がアクセス制限で取得失敗)。
- [事実] QDD 用 1 段遊星の比較研究では、試作品が 6:1 である。Singh, Kapa, Chedda, Kolathaya, "Comparison between External and Internal Single Stage Planetary gearbox actuators for legged robots," arXiv:2506.16356, 2025. https://arxiv.org/abs/2506.16356
- 逆駆動性を左右する要因 [一般的な機械工学の知見。本調査では個別出典なし]
  - 反射慣性: 出力側から見たロータ慣性は減速比の 2 乗倍になる。
  - 摩擦: 出力側から見たモーターのコギングと摩擦は減速比倍になる。
  - 逆効率: 段数を重ねると効率が掛け算で下がる。逆駆動の効率は正駆動より低い。
  - 予圧: バックラッシをなくす予圧 (鋼球式やサイクロイドで使う) は摩擦を増やす。
  - 典型的な QDD の減速比は 6〜10:1 で、1 段で済む範囲である。

## 3. 3D プリント製 QDD・低減速アクチュエータの事例

| 名前 | URL | 方式 | 減速比 | モーター | 外径 | 材料 | 性能数値 |
|---|---|---|---|---|---|---|---|
| Berkeley Humanoid Lite 6512 アクチュエータ | https://arxiv.org/html/2504.17249v1 | 印刷サイクロイド | 本文に記載なし (二次資料では 15:1、未検証) | MAD M6C12 150KV | 記載なし (6811ZZ 軸受で寸法を決定) | PLA | バックラッシ最大 0.0229 rad、減速機効率 約 90%、剛性 319.49 N·m/rad、60 h 耐久試験で摩耗は許容範囲 |
| Berkeley Humanoid Lite 5010 アクチュエータ | 同上 | 印刷サイクロイド | 記載なし | 5010 ドローン用 | 記載なし | FDM (材料の記載なし) | 記載なし |
| TipTapMotor (OpenTorque 系) | https://github.com/DarrenLevine/TipTapMotor | 1 段遊星 | 8:1 | iPower GBM3506H-130T (ジンバル用) | 47 × 47 × 47 mm | 印刷 (材料の記載なし) | 記載なし |
| (8:1) 3D Printed Planetary Gear Actuator Module (Dongil Choi) | https://hackaday.io/project/167656 | 1 段遊星、モジュール 0.65 (歯数 10/30/70) | 8:1 | Tarot 4008 330KV | 遊星側の外径 56 mm | 印刷 (材料の記載なし) | 「バックラッシ極小、逆駆動性良好」(定性のみ)、繰返し精度 SD 0.0058 mm @150 mm |
| Caden Kraft 印刷 QDD | https://www.hackster.io/news/caden-kraft-s-3d-printable-quasi-direct-drive-actuators-cost-as-little-as-40-each-b1f5fa7376d7 | ステータ内蔵の遊星。丸歯形 (pygeartrain で生成) | 記載なし | 自作 36 スロット、印刷ハルバッハ配列 | 記載なし | PA6-GF | 保持トルク 29.4 N·m (コントローラの 50 A 制限で停止) |
| SpryDrive | https://github.com/bbokser/sprydrive | 1 段遊星 (鋼製歯車。ハウジングのみ印刷) | 7:1 | Herlea X8318S / Turnigy 9235 | 記載なし | 鋼歯車、PLA/PETG/PA ハウジング | バックラッシ 約 0.13° |
| Cycloidal QDD (Seojin Jeong) | https://jeongseojin.github.io/blog/2025/qdd-actuator/ | 内蔵サイクロイド、2 枚 180° 位相、ローラ出力 | 10:1 | 自作 8110 ステータ 36N42P | 記載なし | PA-CF12 ハウジング、試作の歯車は CNC アルミ | 最大 8.8 N·m、逆駆動可能 (定性のみ) |
| Modular 3D-printed proprioceptive actuator (Seojin Jeong) | https://jeongseojin.github.io/blog/2025/modular-3d-printed-proprioceptive-actuators | 鋼球要素の波動型減速機 (Mishin Machine の設計を参考) | 記載なし | X4108S / 5010 / ML5210 | 記載なし | 印刷 (材料の記載なし) | 逆駆動可能 (定性のみ) |
| Yoshida ら | DOI 10.1109/SII55687.2023.10039274 | 印刷遊星 | 8:1 (検索結果の抜粋による) | 記載なし | 記載なし | FFF 樹脂 | 静的効率、剛性、耐久を測定 (数値は本文未取得) |
| Bezzini ら | https://www.synapsesocial.com/papers/69df2bece4eeef8a2a6b0c9d | 印刷複合遊星。印刷軌道面と鋼ころ | 30:1 | 記載なし | 記載なし | 印刷樹脂 | 逆駆動トルク 0.42 ± 0.06 N·m、軸受版比で質量 −32 % |

書誌:
- Chi, Liao, Long, Huang, Shao, Nikolic, Li, Sreenath, "Demonstrating Berkeley Humanoid Lite: An Open-source, Accessible, and Customizable 3D-printed Humanoid Robot," RSS 2025, arXiv:2504.17249. 全部品が 200 × 200 × 200 mm に収まる。
- Yoshida, Endo, Okubo, Nabae, "Experimental Evaluation of a Quasi-direct-drive Actuator with a 3D-printed Planetary Gear Reducer," IEEE/SICE SII 2023. DOI 10.1109/SII55687.2023.10039274. 本文は未読 (T2R2 が停止)。
- Okubo, Nabae, Endo, 印刷トロコイド減速機の耐久評価, ROBOMEC 2021 (https://www.jstage.jst.go.jp/article/jsmermd/2021/0/2021_2P3-A04/_article/-char/en)。検索結果の抜粋によれば約 929 時間で破損した。本文は未確認。
- Bezzini, Bassani, Avizzano, 3D 印刷複合遊星減速機, Journal of Mechanical Design。書誌の詳細は未確認。

## 4. 方式比較

評価は主に [推測]。出典があるものは本文中に示す。

| 方式 | 逆駆動性 | バックラッシ | FDM 適性 | 小型化 |
|---|---|---|---|---|
| 1 段遊星 (平歯) | 高い。6〜10:1 で 1 段 | 中程度 | 中程度。小さい歯に応力が集中し、衝撃で折れやすい [事実: Jeong の記事] | 高い。ロータ内に収められる |
| 1 段遊星 (やまば歯) | 高い | 平歯より小さい | 良い。積層方向に強く、噛み合いが滑らか [事実: TNU 研究] | 同上。軸方向の寸法が増える |
| サイクロイド (ピン/ローラ) | 中〜高。10:1 程度なら可。ピンにローラを付けると摩擦が下がる | 小さい (多点で接触) | 高い。荷重が複数の歯に分散し、印刷精度の粗さを許容する [事実: Berkeley Lite] | 中程度。偏心軸受が必要 |
| 鋼球式 (ball reducer) | 予圧の量で決まる。QDD での定量報告は見つからなかった | 予圧すればほぼゼロ [事実: メーカー資料、Terada 1988] | 中程度。溝の精度が効く | 高い。薄型にできる |
| ストレインウェーブ | 低い。高減速比 (30:1 以上) で QDD に向かない | ほぼゼロ | 低い。柔軟な歯車の寿命が問題 | 高い |
| 複合遊星 (Wolfrom) | 低〜中。高減速比で効率が低い。印刷版の逆駆動トルクは 0.42 N·m [事実: Bezzini] | 中程度 | 中程度 | 高い |
| キャプスタン/ベルト | 最も高い。バックラッシなし | ほぼゼロ | 高い (プーリを印刷) | 低い。大径になる |

鋼球式について:
- [事実] 原理は Terada らの「サイクロイドボール減速機」で、2 枚の溝付き円板の間に鋼球を挟む。精密工学会誌 54(11):2101–2106, 1988 (https://www.jstage.jst.go.jp/article/jjspe1986/54/11/54_11_2101/_article/-char/ja)。効率の理論計算法は同誌 63(6):834–838, 1997 (https://www.jstage.jst.go.jp/article/jjspe1986/63/6/63_6_834/_article/-char/ja)。
- [事実] QDD に使って逆駆動性や効率を定量的に報告した事例は、今回の検索では見つからなかった。最も近いのは、鋼球要素の波動型減速機を使った Jeong の試作 (定性のみ) である。

## 5. FDM で歯車を作るときの知見

- [事実] モジュール 0.65 の 1 段遊星 (8:1) を印刷した実例がある (Dongil Choi, https://hackaday.io/project/167656)。0.4 mm ノズルでは、これが実用下限の目安になる。
- [事実] 歯幅の経験則は、モジュールの 12〜20 倍。インフィルは 40〜60 %、材料は PETG か PA を推奨する (MakerWorld のパラメトリックモデル, https://www.makerworld.com/en/models/2847786)。試験データではなく経験則。
- [事実] How to Mechatronics の破壊試験では、インフィル 35 % 以上、壁 5 周以上を推奨している。材料の吸湿や経年でも強度が変わる (https://howtomechatronics.com/how-it-works/how-to-3d-print-gears-the-ultimate-guide/)。
- [事実] PLA, ABS, PETG, ASA を比較した研究では、PLA と ASA の破断荷重が高かった。やまば歯は、滑らかさと耐荷重の釣り合いが最も良かった (https://jst.tnu.edu.vn/jst/article/download/12714/pdf)。
- [事実] PA6-GF のガラス繊維は相手歯車を摩耗させうるため、歯車には素のナイロンがよいという実践者の指摘がある (Caden Kraft, 上記 hackster)。
- [事実] やまば歯は軸方向の推力を打ち消す。ただし軸受は軸方向の荷重も受けられるものにする (https://hackaday.com/2024/06/12/studying-the-finer-points-of-3d-printed-gears)。
- [推測] 歯は寝かせて (軸を Z 方向にして) 印刷すると、歯元が積層方向で割れにくい。

## 6. 推奨 (GM5208-24 と Creality Hi の組み合わせ)

[推測。数値は上の出典値からの概算] 推す方式は、モーターと同軸に重ねる 1 段遊星 (やまば歯) で、減速比は 6〜8:1、外径は 75〜85 mm である。理由は次のとおり。第一に、GM5208-24 の負荷トルクは 0.18〜0.25 N·m と小さく、13.7 Ω の高抵抗巻線なので電流も約 1.75 A @24 V に制限される。したがって歯にかかる荷重は小さく、Mini Cheetah や Berkeley Lite でサイクロイドを選ぶ理由になった衝撃荷重と歯の強度の問題は、このモーターでは起きにくい。第二に、1 段遊星なら逆効率が高い。予圧が要らないので、低トルクモーターでも出力側から回せる。第三に、ジンバル用モーターを 8:1 の印刷遊星と組んだ前例として TipTapMotor (GBM3506H) と Choi (Tarot 4008, モジュール 0.65) がある。出力トルクは、連続 0.2 N·m × 8 × 効率 0.85 ≒ 1.4 N·m、短時間 0.35 N·m × 8 × 0.85 ≒ 2.4 N·m が目安になる。最大トルクはメーカー非公表なので、実測で確かめる必要がある。歯は PETG で、モジュール 0.8〜1.0、歯幅 10〜15 mm とし、0.4 mm ノズルで印刷する。目安の歯数は太陽 12、遊星 24、内歯車 60 で 6:1 となる。軸受と出力キャリアは市販品を使う。モーター外径 59.5 mm と内歯車の歯底、ハウジング肉厚を合わせると、外径はおよそ 80 mm に収まる。造形範囲 260 mm 角に対して余裕は十分ある。サイクロイドは、10:1 以上が必要になったときや衝撃荷重を想定するときの第二候補である。現行の鋼球式は、QDD での効率や逆駆動性の公開データがない。予圧による摩擦が小さいトルクのモーターに効くので、そのまま流用するなら、逆駆動トルクの実測を先に行うべきである。取付穴の寸法は公式に記載がないので、実物で測る必要がある。
