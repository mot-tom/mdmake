A.ヘルプを見る
$ mdmake -h
　　→ヘルプを参照できる

-------------------------------------------------------------------------------------------------------

B.オプションの説明
  (構造条件)
  -c:組成(分子データ作成 & 酸化データ作成)
     ファイル名(convert & stack & cluster)

  -r:組成比(分子データ作成:比率が一定でない、2種以上の混合物に適用できる)

  -s:応力[GPa](分子データ作成:応力値で結晶を伸び縮みさせる際に適用。xx,yy,zz,yz,zx,xyの6値にそれぞれ入力する
                               ※注意:弾性定数が決められていない組成には適用できない)

  -e:歪み[%](分子データ作成:比率で結晶を伸び縮みさせる際に適用。xx,yy,zz,yz,zx,xyの6値にそれぞれ入力する
                               ※応力とは異なり、全ての構造に適用して伸び縮みさせることが可能)

  -n:大きさ(分子データ作成 & 酸化データ作成)
     倍化(convert)
     計算回数(cluster [例]-n "1 10 10" とすると、10回の入れ替え計算を10回繰り返す)

   >:出力(分子データ作成 & convert & stack)


  (周期的境界条件 ※clusterを除く)
  -l:PBC(periodic boundary condition) layer
     指定した方向(x or y or z)の周期的境界条件を無くして層状構造にする
  -w:PBC wired
     指定した方向(x or y or z)以外の周期的境界条件を無くしてワイヤ型構造にする
  -f:PBC free
     全ての方向の周期的境界条件を無くす


  (形状指定)
  -t:形状指定(分子データ作成)
    cuboid        :直方体(デフォルト)
    sphere        :球体
    cylinder      :x方向に延びる円柱(yやzでも境界で接円)
    cylinder_y    :y方向に延びる円柱(xやzでも境界で接円)
    cylinder_z    :z方向に延びる円柱(xやyでも境界で接円)

    wavedcylinder2:x方向に延びる、くびれ(正弦波)が2つある円柱(yやzでも境界で接円)
    wavedcylinder4:x方向に延びる、くびれ(正弦波)が4つある円柱(yやzでも境界で接円)
    sawedcylinder2:x方向に延びる、くびれ(三角波)が2つある円柱(yやzでも境界で接円)
    sawedcylinder4:x方向に延びる、くびれ(三角波)が4つある円柱(yやzでも境界で接円)
    wavedprism2   :x方向に延びる、くびれ(正弦波)が2つある四角柱
    wavedprism4   :x方向に延びる、くびれ(正弦波)が4つある四角柱
    sawedprism2   :x方向に延びる、くびれ(三角波)が2つある四角柱
    sawedprism4   :x方向に延びる、くびれ(三角波)が4つある四角柱

    sinwinding1   :x方向に延びる、うねり(正弦波)が1回あるワイヤー
    sinwinding2   :x方向に延びる、うねり(正弦波)が2回あるワイヤー
    triwinding1   :x方向に延びる、うねり(三角波)が1回あるワイヤー
    triwinding2   :x方向に延びる、うねり(三角波)が2回あるワイヤー
    sqrwinding1   :x方向に延びる、うねり(方形波)が1回あるワイヤー
    sqrwinding2   :x方向に延びる、うねり(方形波)が2回あるワイヤー

　　※断面積が等しくなる構造を生成。また、yとzも同様
    cuboid_xn     :x方向に伸びる円柱(xのみ周期的)
    cylinder_xn   :x方向に延びる四角柱(xのみ周期的)
    prism_xn      :x方向に延びる二等辺三角柱(xのみ周期的)


  (その他)
  -p     :lmpファイルに出力する際にポテンシャルを指定する際に使用

  --seed :seed値指定

  --omp  :OpenMPに使用するプロセッサ数指定

-------------------------------------------------------------------------------------------------------

C.mdmakeのコンパイル方法
  ①Makefile内にあるOBJSに準じた.cppファイルが全てあるか確認する。
　②$ make
　　→エラーが出なければコンパイル成功
　　　(時刻の同期ずれが起きた場合、)

-------------------------------------------------------------------------------------------------------

1.結晶の分子データを作成する
$ mdmake -c Si -n "4 4 4" > Si.mdl
　　→Siの結晶をxyzの3方向に4つ分ずつ作成し、Si.mdlに出力する

組成表("-c "のあとに適用できる組成)
  ※使用可能なオプション(-r,-s)をプログラム名に付記
  [diamond.cpp] -r(単結晶以外),-s
  Si,Ge,C,Sn      :ダイヤモンド構造の単結晶
  SiGe,SiC        :ダイヤモンド構造でランダムに配置された2種混合結晶
  GeSiSn          :ダイヤモンド構造でランダムに配置された3種混合結晶

  [wurtzite.cpp] -r(SiGe,GeSiSnのみ),-s
  2H-Si           :種類の異なる六方晶構造が2層ずつ積み重なって出来た構造(パターンAB)
    (SiC,AlN,BN,BeO)
  3H-SiC          :パターンABC
  4H-Si(SiC)      :パターンABCB
  6H-SiC          :パターンABCACB
  8H-SiC          :パターンABCABACB
  10H-SiC         :パターンABCABCBACB
  12H-SiC         :パターンABCABCACBACB
  15R-SiC         :パターンABCACBCABACABCB
  18H-SiC         :パターンABCACBABCBABCBABCB
  3H-Si           :Siの結晶方向違い,x:(1-10) y:(11-2) z:(111)(=六方晶が3層重なったパターンABCと同一)
  3H-SiGe-alloy   :3H-SiのSiGe版
  3H-GeSiSn-alloy :3H-SiのGeSiSn版

  [tetragonaldia.cpp] -r(Si以外),-s
  3T-Si           :Siの結晶方向違い,x:(110) y:(1-10) z:(001)
  3T-SiGe-alloy   :3T-SiのSiGe版
  3T-GeSiSn-alloy :3T-SiのGeSiSn版

  [orthodia.cpp] -r(Si以外),-s
  3HO-Si          :Siの結晶方向違い,x:(111) y:(11-2) z:(1-10)
  3HO-SiGe-alloy  :3HO-SiのSiGe版
  3HO-GeSiSn-alloy:3HO-SiのSiGe版

  [orthodia2.cpp] -r(Si以外),-s
  3TO-Si          :Siの結晶方向違い,x:(100) y:(01-1) z:(011)
  3TO-SiGe-alloy  :3HO-SiのSiGe版
  3TO-GeSiSn-alloy:3HO-SiのSiGe版

  [confinementdia.cpp]
  Si-in-Ge        :Geに閉じ込められたSiのダイヤモンド構造
  Ge-in-Si        :Siに閉じ込められたGeのダイヤモンド構造

  [superlattice.cpp]
  sl1-SiGe        :SiとGeを交互に1層ごと積層した構造
  sl2-SiGe        :SiとGeを交互に2層ごと積層した構造
  slr-SiGe        :SiとGeをランダムに積層した構造
  
  [diamond_oxide.cpp]
  ac-SiO2(GeO2)   :alpha(低温型)クリストバライト構造
  bc-SiO2(GeO2)   :それぞれ結合角の異なる、Beta(高温型)クリストバライト構造
  bc2-SiO2(GeO2)  
  bc3-SiO2(GeO2)  

  [quartz.cpp]
  aq-SiO2(GeO2)   :alphaクォーツ(石英)構造
  bq-SiO2(GeO2)   :betaクォーツ構造

  [tridymite.cpp]
  at-SiO2(GeO2)   :alphaトリディマイト構造
  bt-SiO2(GeO2)   :betaトリディマイト構造

  [rutile.cpp]
  st-SiO2(GeO2)   :スティショバイト構造
  ru-SnO2(TiO2)   :ルチル構造

  [diamond_comp.cpp]
  g-Si3N4         :γ-窒化シリコン
  Y2O3            :イットリア
  C-ZrO2(SnO2)    :cubic-ジルコニア(酸化スズ)

  [tetra_comp.cpp]
  an-TiO2         :アナターゼ型チタニア
  T-ZrO2          :tetragonal-ジルコニア
  b-Sn            :β-スズ
  
  [rocksalt.cpp] -r(GeSbTeのGe:Sb:Vの比率のみ)
  NaCl            :岩塩(塩化ナトリウム)
  MgO,SrO,TiN     :岩塩型構造
  GeSbTe          :Aサイト=(Ge:Sb:Vacancy=2:2:1)、Bサイト=Teとなる、Ge2Sb2Te5の構造

  [corundom.cpp]
  a-Al2O3         :コランダム型構造のアルミナ

  [dichalcogenide.cpp] -r,-s
  MoS2            :ジカルコゲナイド合金の構造
  MoSX,MoXS2      :一部置換用

  [zincblende.cpp]
  GaAs,InAs       :閃亜鉛鉱型構造
  3C-SiC(SiGe,BN) :種類の異なる立方晶構造が3層ずつ積み重なって出来た構造
  3C-SiCX         :一部置換用

  [fcc.cpp]
  Ni,Au           :面心立方格子型構造
   
-------------------------------------------------------------------------------------------------------

2.アモルファスの分子データを作る(am-)

$ mdmake -c am-Si -n "4 4 4" > am-Si.mdl
　　→Siの結晶をxyzの3方向に4つ分ずつ作成した後にアモルファス化し、am-Si.mdlに出力する

組成表("-c am-"のあとに適用できる組成)
  [oxi_amo.cpp]
  Si,aq-SiO2,bq-SiO2,at-SiO2,bt-SiO2,st-siO2,ac-SiO2,bc-SiO2,
  Ge,bc-GeO2,
  ru-TiO2,a-Al2O3,MgO,SrO

-------------------------------------------------------------------------------------------------------

3.z軸レイヤー情報付きの結晶分子データを作る(zl-,zl1-)

$ mdmake -c zl-Si -n "4 4 4" > zl-Si.mdl
　　→Siの結晶面をxy方向に4格子分作成し、中心からz軸の正負両方向に4層づつ積み上げて、層番号を付加し、zl-Si.mdlに出力する

$ mdmake -c zl1-Si -n "4 4 4" > zl1-Si.mdl
　　→Siの結晶面をxy方向に4格子分作成し、中心からz軸の正方向にのみ4層積み上げて、層番号を付加し、zl1-Si.mdlに出力する

組成表("-c zl-"または"-c zl1-"のあとに適用できる組成)
  [diamond_layer.cpp]Si,
                     Ge(zl-Geのみ)
  [wurtz_layer.cpp]  3H-Si
  [tetla_layer.cpp]  3T-Si
  [ortho_layer.cpp]  3HO-Si
  [ortho2_layer.cpp] 3TO-Si

-------------------------------------------------------------------------------------------------------

4.z軸レイヤー情報付き結晶を、層ごとに酸化する

$ mdmake -c zl-SiO2 -n "4 4 4"
　　→zl-Siで出力したファイルを、z軸両方向にあるそれぞれの表面から層ごとに酸化する。

組成表("-c "のあとに適用できる組成)
  [oxi_layer.cpp]
  zl-SiO2,zl1-SiO2,zl-GeO2
  zl-3H-SiO2,zl1-3H-SiO2
  zl-3T-SiO2,zl1-3T-SiO2
  zl-3HO-SiO2,zl1-3HO-SiO2
  zl-3TO-SiO2,zl1-3TO-SiO2

-------------------------------------------------------------------------------------------------------  

5.ファイルの形式を変更する

$ mdmake --convert lmp -n "2 1 1" -c test.mdl > test.lmp
　　→mdl形式であるtest.mdlを、x方向に2倍の大きさにしてlmp形式(lammps計算用)に変換する。

$ mdmake --convert mdl --lammps SiO2 -c test.final > test.mdl
　　→lmpの計算結果であるtest.final(組成:SiO2)を、mdl形式(lammps計算用)に変換する。

変換形式表("--convert "のあとに適用できる形式)
  [convert.cpp]
  mdl         :final形式やmdl形式のデータをmdl形式に変換する
  mdlZtoX     :mdl形式のデータをy軸で回転させる
  mdlYtoX     :mdl形式のデータをz軸で回転させる
  mdlCUT      :mdl形式のデータで、z軸方向ではみ出た原子を消去する
  mdlFILM     :mdl形式のデータで、z軸の両方向に元の分と同じだけの何もないスペースを確保する(zの大きさは3倍)
  mdlWIRE     :mdl形式のデータで、全軸両方向に元の分と同じだけの何もないスペースを確保する(x,y,zの大きさは全て3倍)
  mdlFIX      :mdl形式のデータで、z=0.00の原子を固定する
  mdlXO2      :mdl形式のデータで、SiO2をSXO2(SXは置換用)に変換する。
  mdlCLEAN    :mdl形式のデータで、最も大きい塊以外の小片をすべて消す
  mdlZRES     :mdl形式のデータで、z軸の両方向に"本来の大きさ"と同じだけの何もないスペースを確保する
  mdlGROUP    :SiGeを含むmdl形式のデータで、隣接する同原子らの塊をグループ分けする。
  mdlINTERFACE:酸化物による界面を含むmdl形式のデータで、界面に存在する原子を減らす。
  lmp     :final形式やmdl形式のデータをlmp形式に変換する
  lmpCIM  :final形式やmdl形式のデータをlmp形式に変換する(CIMポテンシャル用) ※oxi_amoで対象を比較すると分かりやすい
  lmpCMAS :final形式やmdl形式のデータをlmp形式に変換する(CMASポテンシャル用)
  lmpAr   :lammpsでArイオンを照射する際に変換する特殊なlmp形式
  lmpLAY  :z軸レイヤー情報を保存したままlmp形式に変換する
  xyz     :xyz形式に変換する ※VMDやVESTAといったソフトで視覚的に分子を観察できる。
  vasp    :VASP用にファイルを変換する
  alm     :ALAMODE用にファイルを変換する


-------------------------------------------------------------------------------------------------------

6.結晶をくっつける[stack.cpp]

$ mdmake --stack test.mdl -c test2.mdl > test3.mdl
　　→test.mdlとtest2.mdlをz方向でくっつけてtest3.mdlに出力する


-------------------------------------------------------------------------------------------------------

7.クラスターを作る(SiGe限定)

$ mdmake --cluster pxyz -n "1 10 10"-c test.mdl
　　→周期的境界条件を全方向に適用したうえで、10回の入れ替え計算を10回繰り返す
　　(-nの最初の数字"1"は無意味)

変換形式表("--cluster "のあとに適用できる形式)
  [cluster.cpp]
  none  : 周期的境界条件
  px    : x方向に周期的境界条件を設定する
  py    : y方向に周期的境界条件を設定する
  pz    : z方向に周期的境界条件を設定する
  pxy   : x方向とy方向に周期的境界条件を設定する
  pyz   : y方向とz方向に周期的境界条件を設定する
  pxz   : x方向とz方向に周期的境界条件を設定する
  pxyz  : 全方向に周期的境界条件を設定する

  
-------------------------------------------------------------------------------------------------------

8.動径分布関数を計測する

$ mdmake --rdf rdf -n "1 1 10"-c test.mdl > test.rdf
　　→test.mdlを入力にとり、全種類のうち1種類目の原子(convert準拠)を中心に、
　　　同じ1番種類目の原子を対象とする動径分布関数を距離10Å以内で計算する。
　　　(SiO2の場合、Siが1種類目でOが2種類目であり、上記ではSi-SiのRDFを計算する)

$ mdmake --rdf num -n "1 2 2.7"-c test.mdl > test.rdf
　　→test.mdlを入力にとり、1番目の原子(convert準拠)を中心に、
　　　2種類目の原子を対象とする配位数を距離2.7Å以内で計算する。

変換形式表("--cluster "のあとに適用できる形式)
  [cluster.cpp]
  rdf : 動径分布関数を計算する
  num : 配位数を計算する


