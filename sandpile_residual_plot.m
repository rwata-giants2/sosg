function sandpile_residual_plot()
    % 各格子サイズLについてシミュレーション
    L_values = [20, 35, 50];       % 比較したい格子サイズ
    iteration = 15000;              % 時間ステップ
    threshold = 4;            %雪崩が発生する閾値
    %1つのセルの砂粒数がthresholdを超えると雪崩

    % 結果格納用
    totalSand = zeros(length(L_values), iteration); %行列の初期化

    % 各サイズLについて繰り返し
    for idx = 1:length(L_values) %各格子サイズごとに繰り返す
        L = L_values(idx); 
        z = zeros(L, L);           %L×L格子初期状態

        for t = 1:iteration % 1粒砂をランダム位置に追加
            x = randi(L);   %格子のx方向の位置
            y = randi(L);   %格子のy方向の位置
            z(x, y) = z(x, y) + 1; %ランダムに選ばれたセルに砂粒を１つ追加

            %崩壊する時の処理
            unstable = true; %セルに閾値以上の砂粒
            while unstable %閾値以上の砂粒がある限り，崩壊処理を繰り返す
                unstable = false; %セルに閾値未満の砂粒
                [I, J] = find(z >= threshold); %閾値以上のセルを探す
                for k = 1:length(I) %length(I)回ループ、不安定セルを1つずつ処理
                    i = I(k); %k番目の不安定セルの位置iを特定
                    j = J(k); %k番目の不安定セルの位置jを特定
                    z(i,j) = z(i,j) - 4; %セルが崩壊し,砂を4方向にばらまく
                    if i > 1
                        z(i-1,j) = z(i-1,j) + 1; %セルの上方向に砂を1粒渡す
                    end
                    if i < L
                        z(i+1,j) = z(i+1,j) + 1; %セルの下方向に砂を1粒渡す
                    end
                    if j > 1
                        z(i,j-1) = z(i,j-1) + 1; %セルの左方向に砂を1粒渡す
                    end
                    if j < L
                        z(i,j+1) = z(i,j+1) + 1; %セルの右方向に砂を1粒渡す
                    end
                    unstable = true; %砂粒を渡し，崩壊が起こった
                end
            end

            % 現在の砂粒合計を記録
            totalSand(idx, t) = sum(z(:));
        end
    end

    % 結果プロット
    figure;
    hold on;
    colors = {'b', 'g', 'r'};
    for idx = 1:length(L_values)
        plot(1:iteration, totalSand(idx, :), ...
            'Color', colors{idx}, ...
            'DisplayName', sprintf('L=%d', L_values(idx)));
    end
    xlabel('t');
    ylabel('The number of sand grains');
    legend show;
    title('Residual Sand Grains over Time for Various L');
    grid on;
    exportgraphics(gcf, 'Time series of number of sand grains remaining on the lattice with different L in BTW sandpile model.pdf');
end


