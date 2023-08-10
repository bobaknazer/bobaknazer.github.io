function [sumcapacity,rates_SIC,dominant_rates_para,optimal_rates_succ,dominant_matrices_para,optimal_matrices_succ] = compmacweb(H,powervec,figstart)

%Function to generate and plot multiple-access rate tuples for parallel and successive %compute-and-forward. Code is currently limited to L = 2 users and real-valued channels.
%Checks all viable matrices A, which may be quite slow. To restrict the search space,
%use the parameter below to limit the absolute value of the integer entries of possible matrices.

maxint = 5;

[N L] = size(H);
P = diag(powervec);
P1 = P(1,1);
P2 = P(2,2);

F = [chol(inv((inv(P) + H'*H)))]; %Calculate "channel" lattice
                                   %generator F

sumcapacity = 1/2 * log2(det(eye(N) + H*P*H')); %Sum capacity

%SIC corner points
Hcol1 = H(:,1);
Hcol2 = H(:,2);
R_SIC_1a = 1/2 * log2(1 + P1*Hcol1'*inv(eye(N) + P2*Hcol2*Hcol2')*Hcol1);
R_SIC_2a = 1/2 * log2(1 + P2*Hcol2'*Hcol2);
R_SIC_1b = 1/2 * log2(1 + P1*Hcol1'*Hcol1);
R_SIC_2b = 1/2 * log2(1 + P2*Hcol2'*inv(eye(N) + P1*Hcol1*Hcol1')*Hcol2);
rates_SIC = [R_SIC_1a R_SIC_2a; R_SIC_1b R_SIC_2b];

%Capacity region line
capacity_region = [0 R_SIC_2a; rates_SIC; R_SIC_1b 0];

%Find maximum eigenvalue
eigenvalues = eig(eye(2) + P*H'*H);
maxeig = max(eigenvalues);
maxeig = min(maxint^2,maxeig); %truncates search

possible_integers = min(floor(-sqrt(maxeig)),-1):max(ceil(sqrt(maxeig)),1);

[a11ndgrid a12ndgrid a21ndgrid a22ndgrid] = ndgrid(possible_integers);

num_possible_matrices = (length(possible_integers))^4;

a11temp = reshape(a11ndgrid,1,num_possible_matrices);
a12temp = reshape(a12ndgrid,1,num_possible_matrices);
a21temp = reshape(a21ndgrid,1,num_possible_matrices);
a22temp = reshape(a22ndgrid,1,num_possible_matrices);

possible_matrices = zeros(2,2,num_possible_matrices);
possible_matrices(1,1,:) = a11temp;
possible_matrices(1,2,:) = a12temp;
possible_matrices(2,1,:) = a21temp;
possible_matrices(2,2,:) = a22temp;

%Find all parallel computation multiple-access rate pairs

%For each viable integer matrix, generate MAC rates for 
%parallel compute-and-forward. Retain indices of full-rank
%matrices as well as noise variances (to determine dominant 
%solutions afterwards).

rates_para = [];
para_indices = [];
vars_para = [];

for i = 1:num_possible_matrices
  A = possible_matrices(:,:,i);
  if (det(A) ~= 0)
    a1col = A(1,:)';
    a2col = A(2,:)';
  
    %parallel compute-and-forward effective noise variances
    var1_para = sum((F*a1col).^2);
    var2_para = sum((F*a2col).^2);
  
    %three possible cases for admissible mappings
    if ((A(1,1) ~= 0) & (A(1,2) ~= 0))
      %permutation (1,2)
      snr1_para_a = max(1,P1/var1_para);
      snr2_para_a = max(1,P2/max(var1_para,var2_para));
      R1_para_a = 1/2*log2(snr1_para_a);
      R2_para_a = 1/2*log2(snr2_para_a);
      rates_para = [rates_para; R1_para_a R2_para_a];
      para_indices = [para_indices; i];
      vars_para = [vars_para; var1_para var2_para];
    
      %permutation (2,1)
      snr1_para_b = max(1,P1/max(var1_para,var2_para));
      snr2_para_b = max(1,P2/var1_para);
      R1_para_b = 1/2*log2(snr1_para_b);
      R2_para_b = 1/2*log2(snr2_para_b);
      rates_para = [rates_para; R1_para_b R2_para_b];
      para_indices = [para_indices; i];
      vars_para = [vars_para; var1_para var2_para];
    end
  
    if ((A(1,1) ~= 0) & (A(1,2) == 0))
      %permutation (1,2)
      snr1_para = max(1,P1/var1_para);
      snr2_para = max(1,P2/var2_para);
      R1_para = 1/2*log2(snr1_para);
      R2_para = 1/2*log2(snr2_para);
      rates_para = [rates_para; R1_para R2_para];
      para_indices = [para_indices; i];
      vars_para = [vars_para; var1_para var2_para];
    end
    
    if ((A(1,1) == 0) & (A(1,2) ~= 0))
      %permutation (2,1)
      snr1_para = max(1,P1/var2_para);
      snr2_para = max(1,P2/var1_para);
      R1_para = 1/2*log2(snr1_para);
      R2_para = 1/2*log2(snr2_para);
      rates_para = [rates_para; R1_para R2_para];
      para_indices = [para_indices; i];
      vars_para = [vars_para; var1_para var2_para];
    end
  end
end

%Find dominant solutions
%Find rate pairs with maximum R1
R1a_dominant = max(rates_para(:,1));
R1dominant_ind = find([rates_para(:,1) == R1a_dominant]);
%amongst those rate pairs, find maximum R2
[R2a_dominant R2a_subind] = max(rates_para(R1dominant_ind,2));

%Find rate pairs with maximum R2
R2b_dominant = max(rates_para(:,2));
R2dominant_ind = find([rates_para(:,2) == R2b_dominant]);
%amongst those rates pairs, find maximum R1
[R1b_dominant R1b_subind] = max(rates_para(R2dominant_ind,1));
dominant_rates_para = [R1a_dominant R2a_dominant; R1b_dominant R2b_dominant];
dominant_matrices_para = zeros(2,2,2);
dominant_matrices_para(:,:,1) = possible_matrices(:,:,R1dominant_ind(R2a_subind)); dominant_matrices_para(:,:,2) = possible_matrices(:,:,R2dominant_ind(R1b_subind));

%Find all successive computation multiple-access rate pairs

unimodular_indices = [];

%Generate list of indices for viable unimodular matrices.
for i = 1:num_possible_matrices
  if (abs(det(possible_matrices(:,:,i))) == 1)
    unimodular_indices = [unimodular_indices i];
  end
end

num_unimodular_matrices = length(unimodular_indices);

%For each viable unimodular matrix, generate MAC rates for 
%successive compute-and-forward. Retain index
%to associated unimodular matrix for each rate pair.

rates_succ = [];
succ_indices = [];

for i = 1:num_unimodular_matrices
  index = unimodular_indices(i);
  A = possible_matrices(:,:,index);
  a1col = A(1,:)';
  a2col = A(2,:)';
  
  %successive compute-and-forward effective noise variances
  var1_succ = sum((F*a1col).^2);
  nullspace = eye(2) - F*a1col*inv(a1col'*F'*F*a1col)*a1col'*F';
  var2_succ = sum((nullspace*F*a2col).^2);
  
  %three possible cases for admissible mappings
  if ((A(1,1) ~= 0) & (A(1,2) ~= 0))
    %permutation (1,2)
    snr1_succ_a = max(1,P1/var1_succ);
    snr2_succ_a = max(1,P2/max(var1_succ,var2_succ));
    R1_succ_a = 1/2*log2(snr1_succ_a);
    R2_succ_a = 1/2*log2(snr2_succ_a);
    rates_succ = [rates_succ; R1_succ_a R2_succ_a];
    succ_indices = [succ_indices; index];
    
    %permutation (2,1)
    snr1_succ_b = max(1,P1/max(var1_succ,var2_succ));
    snr2_succ_b = max(1,P2/var1_succ);
    R1_succ_b = 1/2*log2(snr1_succ_b);
    R2_succ_b = 1/2*log2(snr2_succ_b);
    rates_succ = [rates_succ; R1_succ_b R2_succ_b];
    succ_indices = [succ_indices; index]; 
  end
  
  if ((A(1,1) ~= 0) & (A(1,2) == 0))
    %permutation (1,2)
    snr1_succ = max(1,P1/var1_succ);
    snr2_succ = max(1,P2/var2_succ);
    R1_succ = 1/2*log2(snr1_succ);
    R2_succ = 1/2*log2(snr2_succ);
    rates_succ = [rates_succ; R1_succ R2_succ];
    succ_indices = [succ_indices; index];
  end
  
  if ((A(1,1) == 0) & (A(1,2) ~= 0))
    %permutation (2,1)
    snr1_succ = max(1,P1/var2_succ);
    snr2_succ = max(1,P2/var1_succ);
    R1_succ = 1/2*log2(snr1_succ);
    R2_succ = 1/2*log2(snr2_succ);
    rates_succ = [rates_succ; R1_succ R2_succ];
    succ_indices = [succ_indices; index];
  end
end

%Find sum-rate optimal successive computation rate pairs
sumrates_succ = sum(rates_succ,2);
max_sumrate_succ= max(sumrates_succ);
max_sumrate_succ_indices = find(sumrates_succ > max_sumrate_succ * (1 - 0.0000001));
max_rates_succ = rates_succ(max_sumrate_succ_indices,:);
%remove repeats
[optimal_rates_succ optimal_succ_subind] = unique(max_rates_succ,'rows','stable');
optimal_matrices_succ = possible_matrices(:,:,max_sumrate_succ_indices(optimal_succ_subind));

%Plot dominant parallel computation solutions (of which there are at most 2) and
%sum-rate optimal successive computation solutions.
figure(figstart)
%hold off
plot(capacity_region(:,1),capacity_region(:,2),'k','linewidth',2)
hold on
plot(dominant_rates_para(:,1),dominant_rates_para(:,2),'mx','linewidth',3,'MarkerSize',12)
plot(optimal_rates_succ(:,1),optimal_rates_succ(:,2),'b*','linewidth',2.5,'MarkerSize',12)
plot(rates_SIC(:,1),rates_SIC(:,2),'ro','linewidth',3,'MarkerSize',14)
axis([0 ceil(R_SIC_1b) 0 ceil(R_SIC_2a)])
xlabel('R_1')
ylabel('R_2')
axis square
legend('Capacity Region','Parallel Comp.','Successive Comp.','SIC')
set(findall(gcf,'type','text'),'FontSize',14)
set(gca,'FontSize',14)

%Plot all nonzero parallel and successive computation solutions.
nonzero_para_ind = find(sum(rates_para,2));
unique_rates_para = unique(rates_para(nonzero_para_ind,:),'rows','stable');
nonzero_succ_ind = find(sumrates_succ);
unique_rates_succ = unique(rates_succ(nonzero_succ_ind,:),'rows','stable');

figure(figstart+1)
%hold off
plot(capacity_region(:,1),capacity_region(:,2),'k','linewidth',2)
hold on
plot(unique_rates_para(:,1),unique_rates_para(:,2),'mx','linewidth',3,'MarkerSize',12)
plot(unique_rates_succ(:,1),unique_rates_succ(:,2),'b*','linewidth',2.5,'MarkerSize',12)
plot(rates_SIC(:,1),rates_SIC(:,2),'ro','linewidth',3,'MarkerSize',14)
axis([0 ceil(R_SIC_1b) 0 ceil(R_SIC_2a)])
xlabel('R_1')
ylabel('R_2')
axis square
legend('Capacity Region','Parallel Comp.','Successive Comp.','SIC')
set(findall(gcf,'type','text'),'FontSize',14)
set(gca,'FontSize',14)