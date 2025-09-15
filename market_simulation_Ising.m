function market_simulation_Ising()
    % パラメータの初期化
    iteration = 50000;   % タイムステップ50000
    M = 5;               % 履歴長（historyは0..5^M-1）
    B = 9;               % 最小ロット q_i = floor(w_i/B)
    C = 3;               % 認知閾値（ΔPの量子化に使用）
    pattern = 5^M;

    marketPrice    = 100; %市場価格初期値(任意)
    preMarketPrice = marketPrice; 
    cognitivePrice = 0; %認知的価格P(t)

    % 市場依存 意思決定パラメータ（ノイズ無し）
    delta = 0.05;        % デッドゾーン |x|<delta -> 0
    alpha = 0.90;        % K の慣性0.9
    beta  = -1.0;        % 群衆が当たったら同調0.3
    theta = 0.70;        % 戦略表の重み0.7
    Kclip = 5;           % K の上下限（暴走防止）

    % 状態
    history      = randi([0, pattern-1]);     % %全履歴パターンからランダムに1つを格納
    Dprev        = 0;                         % D(t-1)
    rprev        = 0;                         % r(t-1)
    players      = {};                        % プレイヤ配列

    % 記録
    playerCounts    = zeros(iteration,1);
    returns         = zeros(iteration,1);
    capital_changes = [];

    %% シミュレーションループ
    for t = 1:iteration
        if mod(t,1000)==0
            disp(['Timestep: ', num2str(t), ' / ', num2str(iteration)]);
        end

        N_current = numel(players);
        E    = 0;      % D(t)式の分母(売買の“符号”×“数量”の合計)
        Qcap = 0;      % 合計sum(q) 正規化に分母

        % 1) 意思決定：a^{±}_i(t) = sign( K_i D(t-1) + theta*S_i(H) )
        for i = 1:N_current
            [a_pm1, ~, q_i, players{i}] = decision_with_market( ...
                players{i}, B, history, Dprev, theta, delta);

            % 集計（action_012は互換用。ここでは使わない）
            E    = E    + a_pm1 * q_i;
            Qcap = Qcap + q_i;
        end

        % 2) 市場の向き D(t) と 価格生成
        if Qcap < 1e-12 %？
            D = 0;
        else
            D = E / Qcap;                  % 市場の向きを[-1,1]に
        end

        % 市場価格の変化と定量化
        if N_current > 0
            deltaP = E / N_current;       % ΔP = (1/N) * sum(a^{±} q)
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
        history = mod(history * 5, pattern) + move;

        % 4) 価格とリターン
        marketPrice    = marketPrice + deltaP;
        r              = log(marketPrice) - log(preMarketPrice); %r(t)=Inp(t)-Inp(t-1)
        preMarketPrice = marketPrice; %市場価格の更新

        % 5) プレイヤー状態の更新・退出処理
        if N_current > 0
            keep_flags = true(1, N_current);
            for i = 1:N_current
                players{i} = open_prices(players{i}, cognitivePrice);
                [players{i}, cap_change] = update(players{i}, B, cognitivePrice);

                if cap_change ~= 0
                    capital_changes(end+1) = cap_change; %#ok<AGROW>
                end
                if players{i}.alter
                    keep_flags(i) = false;
                end
            end
            players = players(keep_flags);
        end

        % 6) 新規参入（1人/t）
        players{end+1} = create_player(B, pattern); %#ok<AGROW>

        % 7) K_i の学習更新（残存プレイヤに適用）
        N_current = numel(players);
        for i = 1:N_current
            players{i}.K = players{i}.b + alpha*players{i}.K + beta*(rprev * Dprev);
            players{i}.K = max(min(players{i}.K, Kclip), -Kclip);
        end
        Dprev = D; 
        rprev = r; 

        % 8) 記録
        playerCounts(t) = numel(players);
        returns(t)      = r;
    end

    %% 保存
    writematrix(capital_changes', sprintf('Ising_capital_changes_B%d.csv', B));

    writematrix(playerCounts, sprintf('Ising_player_counts_B%d.csv', B));
    figure; plot(1:iteration, playerCounts);
    xlabel('t'); ylabel('The number of players');
    exportgraphics(gcf, sprintf('Ising_player_counts_B%d.pdf', B));

    writematrix(returns, sprintf('Ising_return_B%d.csv', B));
    figure; plot(returns);
    xlabel('t'); ylabel('r(t)');
    exportgraphics(gcf, sprintf('Ising_return_B%d.pdf', B));
end

%% ===== プレイヤー生成（K,b を追加）=====
function player = create_player(B, pattern)
    player.settle = false;
    player.wealth = 10*B;
    player.quantity = 0;
    player.tradePeriod = 0;
    player.idlePeriod = 0;
    player.idle = 0;
    player.openPosition = 0;
    player.ongoing = 0;
    player.recommends = randi([0,2], 1, pattern); % 0=売 1=待機 2=買
    player.gain = 0;
    player.openPrice = 0;
    player.alter = false;
    player.switched = false;

    % 市場依存パラメータ
    player.K = 0.1;   % 同調強度の初期値
    player.b = 0.0;   % 個人バイアス
end

%% エージェントの意思決定
function [a_pm1, action, q_i, player] = decision_with_market(player, B, history, Dprev, theta, delta)
    % 戦略表の推奨（0/1/2）→ {-1,0,+1}
    select = player.recommends(history+1);
    S_rec  = (select==2) - (select==0);   % -1/0/+1

    % 市場項＋戦略表のみ（ノイズ無し）
    signal = player.K * Dprev + theta * S_rec;

    % デッドゾーン付き符号 → {-1,0,+1}
    a_pm1 = sgn0(signal, delta);

    % SOSG：同方向保有は保留、方向が変わるときは決済
    if player.ongoing > 0
        player.ongoing = player.ongoing + 1;
        if a_pm1 ~= 0 && a_pm1 ~= (player.openPosition-1)
            % 真逆/反転
            player.tradePeriod = player.ongoing - 1;
            player.settle = true;
            player.ongoing = 0;
        elseif a_pm1 ~= 0
            % 同方向
            a_pm1 = 0;
        end
    elseif a_pm1 ~= 0
        % 新規建て
        player.ongoing = 1;
        player.openPosition = (a_pm1<0)*0 + (a_pm1==0)*1 + (a_pm1>0)*2; % 0/1/2に戻す
        player.quantity = floor(player.wealth / B);
        player.idlePeriod = player.idle;
        player.idle = 0;
    else
        player.idle = player.idle + 1;
    end

    % 互換のため 0/1/2 も返す
    if a_pm1 < 0, action = 0; elseif a_pm1 > 0, action = 2; else, action = 1; end
    q_i = max(player.quantity,0);
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

%% ===== 補助 =====
%function mv = quantize_move(deltaP, C)
    % -2..+2 を 0..4 へ写像（既存互換）
    %if     deltaP < -C, mv = 0;
    %elseif deltaP <  0, mv = 1;
    %elseif deltaP >  C, mv = 4;
    %elseif deltaP >  0, mv = 3;
    %else,               mv = 2;
    %end
%end

function y = sgn0(x, delta)
    y = zeros(size(x));
    y(x >=  delta) =  1;
    y(x <= -delta) = -1;
end