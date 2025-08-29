function plot_players_byB()
    % パラメータの初期化
    iteration = 25000; %timestep50000
    %N = 0; %市場参加者数
    M = 5; %履歴長
    B_values = [9,18,27]; %取引最小単位
    C = 3; %認知閾値
    %S = 1; %戦略数
    pattern = 5^M; %全履歴パターン

    colors = {'b','g','r'};
    legends = {};

    allPlayerCounts = zeros(length(B_values), iteration); % 各Bごとの結果記録

    %% シミュレーション開始
    for b_index = 1:length(B_values)
        B = B_values(b_index);
        disp(['シュミレーション開始：B＝', num2str(B)]);
    
        %初期化
        history = randi([0, pattern-1]); %全履歴パターンからランダムに1つを格納
        cognitivePrice = 0; %認知的価格P(t)
        marketPrice = 100; %市場価格の初期値(任意)
        preMarketPrice = marketPrice; %代入
        players = {}; % 空のプレイヤーリスト
        %r = 0; %価格リターン(最初は０)
        playerCounts = zeros(iteration, 1); % プレイヤー数の記録用
        %capital_changes = [];%資本変動量を記録する配列
        %returns = zeros(iteration, 1);% 出力ファイルの初期化

    
        for t = 1:iteration %tを1からinterationまで
            %進行状況可視化
            if mod(t, 500) == 0
                disp(['Timestep: ', num2str(t), ' / ', num2str(iteration)]);
            end

            sell = 0;
            buy = 0;  %売買の合計数を0で初期設定
            N_current = numel(players);


            % 全プレイヤーの注文を取得
            for i = 1:N_current
                [action, players{i}] = decision(players{i}, B, history);
                if action == 0 %売る
                    sell = sell + players{i}.quantity;
                elseif action == 2 %買う
                    buy = buy + players{i}.quantity;
                end
            end

            %deltaP = (buy - sell) / N; %ΔPの式
            N_current = numel(players); %players配列内の要素数=参加者数
            if N_current > 0
                deltaP = (buy - sell) / N_current;
            else
                deltaP = 0;
            end

            %認知閾値とΔPの比較による価格変動
            %0-4で考えるのは後でhistoryを5進数にするため？
            if deltaP < -C
                move = 0;
            elseif deltaP < 0
                move = 1;
            elseif deltaP > C
                move = 4;
            elseif deltaP > 0
                move = 3;
            else
                move = 2;
            end

            cognitivePrice = cognitivePrice + (move - 2); %認知的価格更新
            history = mod(history * 5, pattern) + move; %履歴を5進数？？？
            marketPrice = marketPrice + deltaP; %ΔPを加算
            %r = log(marketPrice) - log(preMarketPrice); %r(t)=Inp(t)-Inp(t-1)
            preMarketPrice = marketPrice; %市場価格の更新

        % プレイヤー状態の更新と退出処理
        %if N_current > 0
            keep_flags = true(1, N_current); % まず全員を残すフラグを立てる
            for i = 1:N_current
                players{i} = open_prices(players{i}, cognitivePrice);
                %update関数から2つの戻り値を受け取る
                [players{i}, cap_change] = update(players{i}, B, cognitivePrice);
      
                %資本変動があった場合(0でない場合)のみ配列に追加
                %if cap_change ~= 0
                    %capital_changes(end+1) = cap_change;
                %end
                
                if players{i}.alter
                    keep_flags(i) = false; % 退出するプレイヤーのフラグを倒す
                end
            end
            players = players(keep_flags); % フラグが立っているプレイヤーだけを残す
        %end

            %new_p = 1;
            %for j = 1:new_p
            players{end+1} = create_player(B, pattern);
            %end

            playerCounts(t) = numel(players); % プレイヤー数を記録
            %playerCounts(b_index,t) = numel(players); % プレイヤー数を記録
            %returns(t) = r; %リターン
        end

        % 結果を保存
        allPlayerCounts(b_index, :) = playerCounts;
        legends{end + 1} = ['B = ', num2str(B)];
    end

        % プレイヤー数のプロット
        figure;
        hold on;
        colors = {'b','g','r'};
        for i = 1:length(B_values)
            plot(1:iteration, allPlayerCounts(i, :), 'Color', colors{i}, 'LineWidth', 1.5,'DisplayName', sprintf('B=%d', B_values(i)));
        end
        plot(1:iteration, allPlayerCounts);
        xlabel('t');
        ylabel('The number of players');
        grid on;
        exportgraphics(gcf, 'Time series of number of players participating in the gamified market with different B in SOSG.pdf');
end


%% プレイヤーの生成(定義的な)
function player = create_player(B, pattern) %作るか
    player.settle = false; %各戦略は今ポジション(売買)を持ってるか
    player.wealth = 10*B;
    player.quantity = 0; %注文数量
    player.tradePeriod = 0; %戦略の取引期間
    player.idlePeriod = 0; %使ってない戦略の非取引期間
    player.idle = 0; %プレイヤーの活動してない期間
    player.openPosition = 0; %？？
    player.ongoing = 0; %戦略に対するポジション保有
    player.recommends = randi([0,2], 1, pattern); %戦略が履歴パターンに推奨する行動
    player.gain = 0; %累積利得
    player.openPrice = 0; %認知的価格？？
    player.alter = false; %退出フラグ=false
    player.switched = false; %switchフラグ= false
end

%% プレイヤーの意思決定
function [action, player] = decision(player, B, history)
    
        select = player.recommends(history+1);

        if player.ongoing > 0 %ポジションを持っているか否か
            player.ongoing = player.ongoing + 1;
            if select ~= 1 && select ~= player.openPosition %待機以外かつポジションと違う行動
                player.tradePeriod = player.ongoing - 1;
                player.settle = true; %戦略はポジションを閉じる
                player.ongoing = 0;
            elseif select ~= 1
                select = 1;
            end
        elseif select ~= 1
            player.ongoing = 1; %今の戦略がbuy,sellを推奨していれば
            player.openPosition = select; %新規取引
            player.quantity = floor(player.wealth / B); %取引量決定
            player.idlePeriod = player.idle; %何もしてない期間を記録更新
            player.idle = 0; %リセット
        else
            player.idle = player.idle + 1; %待機
        end
        
        action = select; %行動を使用中の戦略に従う
end

%% 初期価格記録
function player = open_prices(player, price)
        if player.ongoing == 1
            player.openPrice = price;
        end
end

%% 状態更新
function [player, capital_change] = update(player, B, price)
    capital_change = 0;
    
    if player.settle
            %認知的利得と実際の資産変動額を計算
            gain = (player.openPosition - 1) * (price - player.openPrice);
            delta_w = player.quantity * gain;

            player.gain = player.gain + gain;
            player.wealth = player.wealth + delta_w;

            %Fig.6:資本変動量を計算して戻り値に設定
            capital_change = floor(delta_w / B);
            
            if player.wealth < B
                player.alter = true;
            end
            player.settle = false;
     end
end