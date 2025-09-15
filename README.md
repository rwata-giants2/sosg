# sosg
 Self-organized Speculation Gameの試験的なプログラミングです．

market_simulation.m : 大元となるSOSGの市場シミュレーションのコードです

plot_players_byB : Bの値による市場のスケールの違いを表すコード

sandpile_residual_plot.m : 砂山モデルにおける格子の大きさによる残留砂粒の違いを表すコード

market_sosg_ising_2.m : SOSGをn×n格子にし，Isingモデルの相互作用を導入したもの

market_sosg_ising_fast.m : 隣接プレイヤーの把握方法を変更し、market_sosg_ising_2.mを高速化したもの

market_simulation_Ising.m : SOSGにIsingモデルの相互作用を少し改変し市場の向きや戦略表の遵守具合を導入した
