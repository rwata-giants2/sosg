%%%==========================================================%%%
%%%==========-======market_sosg_ising_fast.m================%%%
%%%=================SOSGにIsing相互作用を導入=================%%%
%%%====隣接するプレイヤーの行動も影響を受ける仕組みを用いる====%%%
%%%==========履歴・戦略表の参照有・外部ニュースの導入有========%%%
%%%==========================================================%%%
%%%==========================================================%%%

function market_sosg_ising_fast(beta)

    if nargin < 1 || isempty(beta)
        beta = 0.50;
    end

    %% パラメータ設定
    %if nargin < 1
        n = 50;
    %end

    iteration = 20000; %timestep50000
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
    capital_changes_prealloc = zeros(iteration * 2, 1); 
    cc_idx = 1; % capital_changesに書き込むためのインデックス
    players = {}; % 空のプレイヤーリスト
    returns = zeros(iteration, 1);% 出力ファイルの初期化

    %Isingモデルの相互作用パラメータ
    delta = 0.05;    % デッドゾーン |x|<delta -> 0
    alpha = 0.7;        % K の慣性
    %beta  = 0.5;        % 群衆が当たったら同調
    theta = 0.5;        % 戦略表の重み
    b_imax = 0.3;
    Kclip = [-3, 3];  % K の上下限（暴走防止）

    lambda = 5;         % Ising価格更新係数 λ
    sigma_max = 0.05;   % 外部ニュース感度 σ_i の上限（例）
    C_V        = 0.01;   % 私的雑音の基準分散（各エージェントに固定）

    %格子
    gridN = n; 
    occupied = false(gridN, gridN);
    player_id_grid = zeros(gridN, gridN); %追加: IDマップ行列を初期化
    D_prev = 0;      %D(t-1)
    r_prev = 0;      %r(t-1)
    G_prev        = 0;                         % G(t-1) の記憶（K学習用）

    
    %% シミュレーションループ
    for t = 1:iteration

        if mod(t, 1000) == 0
            fprintf('t = %d / %d, N = %d\n', t, iteration, numel(players));
        end

        g_t = randn(); % 毎ステップ、新しい外部ニュースを生成
        
        N_current = numel(players);
        playerCounts(t) = N_current;

        %Qcap = 0;      % 合計sum(q) 正規化に分母

        %プレイヤーがいない時は新規参入
        if N_current == 0
            % 新規参入
            [players, occupied, player_id_grid] = enter_market(players, occupied,player_id_grid, B, pattern, b_imax, gridN, C_V, sigma_max);
            %returns(t) = 0;
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
                hloc = (players{i}.K * Es(i)) + theta * s_strat + players{i}.sigma * g_t + players{i}.seps  * randn();

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
            Asum = 0;      % Σ a_i(t) （Ising価格更新

            for i = 1:N_current
                %[action_pm1, players{i}] = decision(players{i}, B, history);
                %if action_pm1 == 2       % buy
                [act012, players{i}] = decision(players{i}, B, history);     % [FIX: 戻り値 0/1/2 に統一]
               
                if act012 == 2
                    buy  = buy  + players{i}.quantity;
                    Qcap = Qcap + players{i}.quantity;
                    Asum = Asum + 1; % 行動が「買い(+1)」なのでAsumに+1

                elseif act012 == 0       % sell
                    sell = sell + players{i}.quantity;
                    Qcap = Qcap + players{i}.quantity;
                    Asum = Asum - 1; % 行動が「売り(-1)」なのでAsumに-1
                end

                % 市場価格の変化と定量化
                if N_current > 0
                    deltaP = (buy - sell) / N_current;
                else
                    deltaP = 0;
                end
    
            end

            % 市場の向き D(t) = (buy - sell)/Σq
            if Qcap > 0
                D = (buy - sell) / Qcap;
            else
                D = 0;
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

            cognitivePrice = cognitivePrice + (move - 2);   %認知的価格更新
            history = mod(history * 5 + move, pattern);
    
            % 4) 価格とリターン
            %marketPrice    = marketPrice + deltaP;
            %marketPrice = max(marketPrice + deltaP, 1e-8);
            %r = log(marketPrice) - log(preMarketPrice);     %r(t)=Inp(t)-Inp(t-1)
            if N_current > 0
                r = (1/(lambda * N_current)) * Asum;      % << NEW: r(t) = (1/(λ N)) Σ a_i(t)
            else
                r = 0;
            end
            marketPrice    = preMarketPrice * exp(r);      % << NEW: p(t) = p(t-1) * exp(r)
            preMarketPrice = marketPrice; %市場価格の更新

            % 5) K更新（b_i + α*K + β*r_prev*G_prev）
            for i = 1:N_current
                Ki = players{i}.b + alpha * players{i}.K + beta * (r_prev * G_prev);
                players{i}.K = clip(Ki, Kclip);
            end
            r_prev = r; 
            G_prev = g_t;
            %D_prev = D;
    
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

                %プレイヤー退出・圧縮後にIDマップを再構築
                %player_id_grid = zeros(gridN, gridN); % いったんリセット
                %for i = 1:numel(players)
                    %pos = players{i}.pos;
                    %player_id_grid(pos(1), pos(2)) = i;
                %end
                [player_id_grid, occupied] = rebuild_maps(players, gridN);
            end
        end
    
            % 7) 新規参入（1人/t）
            %players{end+1} = create_player(B, pattern);
            [players, ~, ~] = enter_market(players, occupied, player_id_grid, B, pattern, b_imax, gridN, C_V, sigma_max);

            [player_id_grid, occupied] = rebuild_maps(players, gridN); %追加

            % 1) occupied と ID マップの一致
            if any(occupied(:) ~= (player_id_grid(:) > 0))
                error('occupied と player_id_grid の不整合を検出');
            end

            % 2) ID マップの範囲チェック
            ids = player_id_grid(player_id_grid>0);
            if ~isempty(ids) && (min(ids)<1 || max(ids)>numel(players))
                error('ID マップが配列添字の範囲外を指しています');
            end

            % 3) 一意性チェック（重複 ID が無いか）
            if numel(unique(ids)) ~= numel(ids)
                error('ID マップに重複 ID が存在します');
            end


    
            % 8) 記録
            playerCounts(t) = numel(players);
            returns(t)      = r;
    end


        %tag = sprintf('n%d_B%d_delta%.1f_theta%.1f_beta%.1f_b_imax%.1f', n, B, delta, theta, beta, b_imax);
        tag = num2str(beta, '%.15g');

        % 資産変動を保存
        %writematrix(capital_changes', sprintf('capital_changes_%s.csv', tag));
        capital_changes = capital_changes_prealloc(1:cc_idx-1);
        writematrix(capital_changes, sprintf('lattice_capital_changes_%s.csv', tag));
        
        % プレイヤー数を保存・描画
        writematrix(playerCounts, sprintf('lattice_player_counts_%s.csv', tag));
        figure;
        plot(1:iteration, playerCounts);
        xlabel('t');
        ylabel('The number of players');
        exportgraphics(gcf, sprintf('lattice_player_counts_%s.pdf', tag));
        
        % リターンを保存・描画
        writematrix(returns, sprintf('lattice_return_%s.csv', tag));
        figure;
        plot(returns);
        xlabel('t');
        ylabel('r(t)');
        exportgraphics(gcf, sprintf('lattice_return_%s.pdf', tag));

end

%% === 意思決定 ===
function [act012, player] = decision(player, B, history)
    select = player.recommends(history+1);  % 0/1/2 = sell/hold/buyに変更
    %action_pm1 = 0; % デフォルトはホールド
    act012 = 1;

    if player.ongoing > 0                                                           %自分のポジションがある
        %player.ongoing は「ポジションを持っている期間（ステップ数）」
        player.ongoing = player.ongoing + 1;                                        %ポジションを保持中なので経過期間を +1
        dir = sign(player.s_prev);                                                  % Ising決定（近傍 + ニュース + 私的ノイズ）
        if dir ~= 0 && dir ~= (player.openPosition - 1)
            player.tradePeriod = player.ongoing - 1;
            player.settle = true;  player.ongoing = 0;
            act012 = 1;  
            return;                              % 決済のみ
        else
            act012 = 1;  return;                 % 方向が同じならホールド
        end

    elseif select ~= 1

                dir = sign(player.s_prev);
                q   = floor(player.wealth / B);
                if dir == 0 || q <= 0
                    act012      = 1;
                    player.idle = player.idle + 1;
                else
                    if dir > 0
                        player.openPosition = 2; act012 = 2;
                    else
                        player.openPosition = 0; act012 = 0;
                    end

                    player.quantity   = q;
                    player.ongoing    = 1;
                    player.idlePeriod = player.idle;
                    player.idle       = 0;
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
    nbrs = [player_pos + [-1 0]; 
            player_pos + [1 0]; 
            player_pos + [0 -1]; 
            player_pos + [0 1]];

    for k=1:4
        x=nbrs(k,1);
        y=nbrs(k,2);

        if x>=1 && x<=n && y>=1 && y<=n
            % IDマップを使って近傍プレイヤーのIDを瞬時に取得
            neighbor_id = player_id_grid(x, y);

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
function player = create_player(B, pattern, position, b_imax, C_V, sigma_max)
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
    player.b = rand() * b_imax;
    player.K = 0.1;                      % 初期 k_i
    player.s_prev = 0;                    % 初期スピン期待
    player.sigma = rand()*sigma_max;     % σ_i ~ U(0, σ_max)
    player.seps   = C_V + 0.1*rand();      % 私的雑音の標準偏差 s_{ε,i}
end

%% === 新規プレイヤー参入 ===
function [players, occupied, player_id_grid] = enter_market(players, occupied, player_id_grid, B, pattern, b_imax, gridN, C_V, sigma_max)
    % 格子の空きを探す
    available_indices = find(~occupied);
    if ~isempty(available_indices)
        idx_to_fill = available_indices(randi(numel(available_indices)));
        [r_idx, c_idx] = ind2sub([gridN, gridN], idx_to_fill);
        
        position = [r_idx, c_idx];

        % <<< 追加: 新しいプレイヤーのID（配列の末尾のインデックス）を取得
        new_player_id = numel(players) + 1;


        players{new_player_id} = create_player(B, pattern, position, b_imax, C_V, sigma_max);
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

function [player_id_grid, occupied] = rebuild_maps(players, gridN)
    player_id_grid = zeros(gridN, gridN);
    occupied       = false(gridN, gridN);
    for i = 1:numel(players)
        pos = players{i}.pos;   % [row, col]
        if ~isempty(pos)
            r = pos(1); c = pos(2);
            if r>=1 && r<=gridN && c>=1 && c<=gridN
                player_id_grid(r, c) = i;   % 「配列添字＝ID」で再割当
                occupied(r, c)       = true;
            end
        end
    end
end