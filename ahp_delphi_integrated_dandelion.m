function ahp_delphi_integrated_dandelion()
    
    % 模拟 Delphi 过程的最终输出
    % 假设 C 矩阵有 3 位专家的判断
    C_expert1 = [1, 0.7; 1.4286, 1]; 
    C_expert2 = [1, 0.9; 1.1111, 1];
    C_expert3 = [1, 0.8; 1.2500, 1];
    C_experts = {C_expert1, C_expert2, C_expert3};

    % 假设 S2 矩阵有 2 位专家的判断
    S2_expert1 = [1, 1.3333; 0.75, 1];
    S2_expert2 = [1, 1.6667; 0.6, 1];
    S2_experts = {S2_expert1, S2_expert2};

    % 假设I22矩阵有4位专家的判断
    I22_expert1 = [1, 2.0; 0.5, 1];
    I22_expert2 = [1, 2.5; 0.4, 1];
    I22_expert3 = [1, 3.0; 0.3333, 1];
    I22_expert4 = [1, 1.8; 0.5556, 1];
    I22_experts = {I22_expert1, I22_expert2, I22_expert3, I22_expert4};

    % 调用子函数geometric_mean_aggregation来聚合Delphi结果
    C = geometric_mean_aggregation(C_experts); 
    S1 = [1, 1; 1, 1]; 
    S2 = geometric_mean_aggregation(S2_experts);
    
    I11 = [1, 1.222; 0.818, 1]; 
    I12 = [1, 1.5; 0.667, 1];   
    I21 = [1, 0.667; 1.5, 1];   
    I22 = geometric_mean_aggregation(I22_experts);

    % AHP计算(调用ahp子函数)
    
    [Cw, Ccr] = ahp(C);
    [S1w, S1cr] = ahp(S1); [S2w, S2cr] = ahp(S2);
    [I11w, I11cr] = ahp(I11); [I12w, I12cr] = ahp(I12);
    [I21w, I21cr] = ahp(I21); [I22w, I22cr] = ahp(I22);
    
    % 综合权重合成
    
    w = [
        Cw(1)*S1w(1)*I11w(1),
        Cw(1)*S1w(1)*I11w(2),
        Cw(1)*S1w(2)*I12w(1),
        Cw(1)*S1w(2)*I12w(2),
        Cw(2)*S2w(1)*I21w(1),
        Cw(2)*S2w(1)*I21w(2),
        Cw(2)*S2w(2)*I22w(1),
        Cw(2)*S2w(2)*I22w(2)
    ];

    
    disp('指标层综合权重：');
    factors = {'土壤改良','传粉服务','资源竞争','化感作用',...
               '食用药用','产业链就业','农业干扰','健康成本'};
    for i=1:8
        fprintf('%d. %-8s 权重：%.4f\n',i,factors{i},w(i));
    end
    fprintf('权重总和：%.4f\n',sum(w));
end

% 子函数

function aggregated_matrix = geometric_mean_aggregation(expert_matrices)
    % 计算多位专家判断矩阵的几何平均值
    
    n = size(expert_matrices{1}, 1);
    num_experts = length(expert_matrices);
    aggregated_matrix = ones(n, n);
    
    for i = 1:n
        for j = 1:n
            product = 1;
            for k = 1:num_experts
                product = product * expert_matrices{k}(i, j);
            end
            aggregated_matrix(i, j) = product^(1/num_experts);
        end
    end
end

function [w, cr] = ahp(mat)
    [V,D] = eig(mat);
    max_eig = max(diag(D));
    w = real(V(:,find(diag(D)==max_eig,1)));
    w = w/sum(w);  
      
    n = size(mat,1);
    if n==1
        cr = 0;
    else
        ci = (max_eig - n)/(n-1);
        ri = [0,0,0.58,0.90,1.12,1.24,1.32,1.41,1.45];
        cr = (n<=length(ri))*(ci/ri(n)) + (n>length(ri))*1;
    end
end