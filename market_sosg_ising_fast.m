function market_sosg_ising_fast(n)

    %% パラメータ設定
    if nargin < 1
        n = 50;
    end

    iteration = 50000; %timestep50000
    M = 5; %履歴長
    B = 9; %取引最小単位
    C = 3; %認知閾値
    pattern = 5^M; %全履歴パターン
    history = randi([0, pattern-1]); %全履歴パターンからランダムに1つを格納
    cognitivePrice = 0; %認知的価格P(t)
    marketPrice = 100; %市場価格の初期値(任意)

    %変数の初期化
    preMarketPrice = marketPrice; %代入
    r = 0; %価格リターン(最初は０)
    playerCounts = zeros(iteration, 1); % プレイヤー数の記録用
    
    %配列の初期化
    %capital_changes = [];%資本変動量を記録する配列
    capital_changes_prealloc = zeros(iteration * 2, 1); % 넉넉한 사이즈를 미리 확보
    cc_idx = 1; % capital_changesに書き込むためのインデックス
    players = {}; % 空のプレイヤーリスト
    returns = zeros(iteration, 1);% 出力ファイルの初期化

    %Isingモデルの相互作用パラメータ
    delta = 0.05;    % デッドゾーン |x|<delta -> 0
    alpha = 0.2;        % K の慣性
    beta  = 1.0;        % 群衆が当たったら同調
    theta = 1.0;        % 戦略表の重み
    b_imax = 1.0;
    Kclip = [-5, 5];  % K の上下限（暴走防止）

    %格子
    gridN = n; 
    occupied = false(gridN, gridN);
    player_id_grid = zeros(gridN, gridN); %追加: IDマップ行列を初期化
    D_prev = 0;      %D(t-1)
    r_prev = 0;      %r(t-1)

    
    %% シミュレーションループ
    for t = 1:iteration

        if mod(t, 1000) == 0
            fprintf('t = %d / %d, N = %d\n', t, iteration, numel(players));
        end
        
        N_current = numel(players);
        playerCounts(t) = N_current;

        %Qcap = 0;      % 合計sum(q) 正規化に分母

        %プレイヤーがいない時は新規参入
        if N_current == 0
            % 新規参入
            [players, occupied, player_id_grid] = enter_market(players, occupied,player_id_grid, B, pattern, b_imax, gridN);
            returns(t) = 0;
            continue;
        end


        if N_current > 0
            %Es = s_j(t-1)
            Es = zeros(N_current,1);

            %高速化のため、事前に全プレイヤーの情報をベクトル化
            all_positions = zeros(N_current, 2); % (N人 x 2次元) の行列を事前確保
            all_s_prev = zeros(1, N_current);    % (1 x N人) のベクトルを事前確保

            for i = 1:N_current
                % player{i}が空でないことを確認してから情報を取得
                if ~isempty(players{i})
                    all_positions(i, :) = players{i}.pos;
                    all_s_prev(i) = players{i}.s_prev;
                end
            end

            for i = 1:N_current
                Es(i) = local_expect_spin_fast(all_positions(i,:), all_s_prev, player_id_grid, gridN);
            end

            %
            for i=1:N_current
                % 戦略表の現時刻推奨をスピンに写像 s_i(H_M(t)) ∈ {-1,0,+1}
                %action_pm1 = players{i}.recommends(history+1);        % 0/1/2
                %s_start = rec2spin(action_pm1); %s_start = s_i(H_M(t))

                act012 = players{i}.recommends(history+1);         % 0/1/2  [FIX: 変数名]
                s_strat = rec2spin(act012);

               
                % sign内=K_i * Es + θ * s_strat
                %hloc = sum( players{i}.K * Es(i) ) + theta * s_start;
                hloc = (players{i}.K * Es(i)) + theta * s_strat;

                % デッドゾーン適用
                if abs(hloc) < delta
                    players{i}.s_prev = 0;
                else
                    players{i}.s_prev = sign(hloc);               % -1/0/+1
                end
                
                % 状態更新
                %players{i}.s_prev = sign(hloc);

                %k_i(t)=
                %k_i(t) = players{i}.b + alpha * players{i}.K +beta*(r_prev*D_prev);
                %players{i}.K = clip(k_i(t), Kclip);
            end

            % 2) 市場の向き D(t) と 価格生成
            buy = 0; 
            sell = 0; 
            Qcap = 0;

            for i = 1:N_current
                %[action_pm1, players{i}] = decision(players{i}, B, history);
                %if action_pm1 == 2       % buy
                [act012, players{i}] = decision(players{i}, B, history);     % [FIX: 戻り値 0/1/2 に統一]
                if act012 == 2
                    buy  = buy  + players{i}.quantity;
                    Qcap = Qcap + players{i}.quantity;
                %elseif action_pm1 == 0   % sell
                elseif act012 == 0       % sell
                    sell = sell + players{i}.quantity;
                    Qcap = Qcap + players{i}.quantity;
                end
            end

            % 市場の向き D(t) = (buy - sell)/Σq
            if Qcap > 0
                D = (buy - sell) / Qcap;
            else
                D = 0;
            end
    
            % 市場価格の変化と定量化
            if N_current > 0
                deltaP = (buy - sell) / N_current;
            else
                deltaP = 0;
            end
    
            % 3) 認知側の更新
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

            cognitivePrice = cognitivePrice + (move - 2);%認知的価格更新
            history = mod(history * 5 + move, pattern);
    
            % 4) 価格とリターン
            %marketPrice    = marketPrice + deltaP;
            marketPrice    = max(marketPrice + deltaP, eps);
            r              = log(marketPrice) - log(preMarketPrice); %r(t)=Inp(t)-Inp(t-1)
            preMarketPrice = marketPrice; %市場価格の更新

            % 5) K更新（b_i + α*K + β*r_prev*D_prev）
            for i = 1:N_current
                Ki = players{i}.b + alpha * players{i}.K + beta * (r_prev * D_prev);
                players{i}.K = clip(Ki, Kclip);
            end
            r_prev = r; 
            D_prev = D;
    
            % 6) プレイヤー状態の更新・退出処理
            if N_current > 0
                keep_flags = true(1, N_current);
                for i = 1:N_current
                    players{i} = open_prices(players{i}, cognitivePrice);
                    [players{i}, cap_change] = update_wealth(players{i}, B, cognitivePrice);
    
                    if cap_change ~= 0
                        %capital_changes(end+1) = cap_change;
                        capital_changes_prealloc(cc_idx) = cap_change;
                        cc_idx = cc_idx + 1;
                    end
                    if players{i}.alter
                        keep_flags(i) = false;
                        pos = players{i}.pos;
                        occupied(pos(1), pos(2)) = false;
                    end
                end
                players = players(keep_flags);

                %追加: プレイヤー退出・圧縮後にIDマップを再構築
                player_id_grid = zeros(gridN, gridN); % いったんリセット
                for i = 1:numel(players)
                    pos = players{i}.pos;
                    player_id_grid(pos(1), pos(2)) = i;
                end
            end
        end
    
            % 7) 新規参入（1人/t）
            %players{end+1} = create_player(B, pattern);
            [players, occupied, player_id_grid] = enter_market(players, occupied, player_id_grid, B, pattern, b_imax, gridN);

    
            % 8) 記録
            playerCounts(t) = numel(players);
            returns(t)      = r;
    end


        tag = sprintf('n%d_B%d_theta%.1f_beta%.1f', n, B, theta, beta);

        % 資産変動を保存
        %writematrix(capital_changes', sprintf('capital_changes_%s.csv', tag));
        capital_changes = capital_changes_prealloc(1:cc_idx-1);
        writematrix(capital_changes, sprintf('capital_changes_%s.csv', tag));
        
        % プレイヤー数を保存・描画
        writematrix(playerCounts, sprintf('player_counts_%s.csv', tag));
        figure;
        plot(1:iteration, playerCounts);
        xlabel('t');
        ylabel('The number of players');
        exportgraphics(gcf, sprintf('player_counts_%s.pdf', tag));
        
        % リターンを保存・描画
        writematrix(returns, sprintf('return_%s.csv', tag));
        figure;
        plot(returns);
        xlabel('t');
        ylabel('r(t)');
        exportgraphics(gcf, sprintf('return_%s.pdf', tag));

end

%% === 意思決定 ===
function [act012, player] = decision(player, B, history)
    select = player.recommends(history+1);  % 0/1/2 = sell/hold/buyに変更
    %action_pm1 = 0; % デフォルトはホールド
    act012 = 1;

    if player.ongoing > 0 %自分のポジションがある
        %player.ongoing は「ポジションを持っている期間（ステップ数）」
        player.ongoing = player.ongoing + 1; %ポジションを保持中なので経過期間を +1
        if select ~= 1 && select ~= player.openPosition %今の推奨が1ではなく，かつ現在の方向とも逆
            % 逆向き推奨が出たらこのステップは決済のみ（新規は出さない）
            player.tradePeriod = player.ongoing - 1;
            player.settle = true; 
            player.ongoing = 0; 
            %action_pm1 = 1; %ホールド
            %action_pm1 = -sign(player.openPosition - 1); % 買いポジション(2)なら-1, 売りポジション(0)なら+1
            act012 = 1;
            return
        end
        % それ以外はポジションをホールド
        act012 = 1;

    elseif select ~= 1
        % 新規建て：方向は s_prev の符号で決める
        %direction = player.s_prev;              % -1/0/+1

        %if direction ~= 0 && floor(player.wealth / B) > 0
            %if direction >= 0
                %player.openPosition = 2; 
            %else 
                %player.openPosition = 0; 
            %end
            %player.quantity = floor(player.wealth / B);
            %player.ongoing = 1;
            %player.idlePeriod = player.idle;
            %player.idle = 0;

            %action_pm1 = direction;

            dir = sign(player.s_prev);                              % -1/0/+1
            q   = floor(player.wealth / B);

            if dir ~= 0 && q > 0
                if dir >= 0
                    player.openPosition = 2; 
                    act012 = 2;
                else
                    player.openPosition = 0; 
                    act012 = 0;
                end
                player.quantity   = q;
                player.ongoing    = 1;
                player.idlePeriod = player.idle; 
                player.idle = 0;
            else
                % 意欲ゼロ or 資金不足なら何もしない
                player.idle = player.idle + 1;
            end
        else
            player.idle = player.idle+1;
    end
end

%% === 近傍期待：上下左右の s_prev の「平均」 ===
function Es = local_expect_spin_fast(player_pos, player_s_prev, player_id_grid, n)
    % player_pos: 対象プレイヤーの座標 [x, y]
    % player_s_prev: 全プレイヤーの前期スピンをまとめたベクトル
    % player_id_grid: 「誰がどこにいるか」のIDマップ行列
    % n: 格子サイズ

    %position = players{i}.pos; 
    Es_sum = 0; 
    cnt = 0;
    nbrs = [player_pos + [ -1 0]; 
            player_pos + [1 0]; 
            player_pos + [0 -1]; 
            player_pos + [0 1]];

    for k=1:4
        x=nbrs(k,1);
        y=nbrs(k,2);

        if x>=1 && x<=n && y>=1 && y<=n
            % IDマップを使って近傍プレイヤーのIDを瞬時に取得
            neighbor_id = player_id_grid(x, y);
            %for j=1:numel(players)
                %if all(players{j}.pos==[x y])
                    %Es_sum = Es_sum + players{j}.s_prev; 
                    %cnt=cnt+1; 
                    %break;


            % IDが0でなければ（＝誰かがいれば）
            if neighbor_id > 0
                Es_sum = Es_sum + player_s_prev(neighbor_id);
                cnt = cnt + 1;
            end
        end
    end
    if cnt > 0
        Es = Es_sum / cnt; % 平均を計算
    else
        Es = 0; % 近傍に誰もいなければ0
    end
end
    %Es = (cnt==0) * 0 + (cnt>0) * (Es_sum/cnt);   % 平均
%end

%% === プレイヤ生成（SOSG+Ising拡張） ===
function player = create_player(B, pattern, position, b_imax)
    player.settle=false; 
    player.wealth=10*B; 
    player.quantity=0;
    player.tradePeriod=0; 
    player.idlePeriod=0; 
    player.idle=0;
    player.openPosition=0; 
    player.ongoing=0;
    player.recommends = randi([0,2],1,pattern);   % SOSGの戦略表
    player.gain=0; 
    player.openPrice=0; 
    player.alter=false; 
    player.switched=false;
    player.pos = position;

    % Ising(改) パラメータ
    player.b = rand()*b_imax;
    player.K = player.b;                      % 初期 k_i
    player.s_prev = 0;                    % 初期スピン期待
end

%% === 新規プレイヤー参入 ===
function [players, occupied, player_id_grid] = enter_market(players, occupied, player_id_grid, B, pattern, b_imax, gridN)
    % 格子の空きを探す
    available_indices = find(~occupied);
    if ~isempty(available_indices)
        idx_to_fill = available_indices(randi(numel(available_indices)));
        [r_idx, c_idx] = ind2sub([gridN, gridN], idx_to_fill);
        
        position = [r_idx, c_idx];

        % <<< 追加: 新しいプレイヤーのID（配列の末尾のインデックス）を取得
        new_player_id = numel(players) + 1;


        players{new_player_id} = create_player(B, pattern, position, b_imax);
        occupied(r_idx, c_idx) = true;

        % 追加: IDマップを更新
        player_id_grid(r_idx, c_idx) = new_player_id;
    end
end

%% === 認知価格のオープン記録（SOSG） ===
function player = open_prices(player, cognitivePrice)
    if player.ongoing==1 
        player.openPrice = cognitivePrice; 
    end
end

%% === 決済・資産更新（SOSG） ===
function [player, cap_change] = update_wealth(player, B, cognitivePrice)
    cap_change=0;
    if player.settle
        gain = (player.openPosition - 1) * (cognitivePrice - player.openPrice);  % 0->-1, 2->+1
        delta_w = player.quantity * gain;
        player.gain = player.gain + gain;
        player.wealth = player.wealth + delta_w;
        cap_change = floor(delta_w/B);
        if player.wealth < B
            player.alter = true; 
        end
        player.settle = false;
    end
end

%% === 補助 ===
function y=clip(x,ab)
    y = max(ab(1), min(ab(2), x)); 
end

function s = rec2spin(rec)   % 0/1/2 -> -1/0/+1
    if rec==0
        s=-1; 
    elseif rec==2
        s=+1; 
    else 
        s=0; 
    end
end