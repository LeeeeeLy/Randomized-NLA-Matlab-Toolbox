clear all
%% set up A
m = 3000; 
n = 300;  
max = 25;         
r = 15;          
p = 5;             
gap = 2;
  
A = controlledgap(m,n,r,gap);   
[U,S,V] = svds(A,max+1);
s = diag(S);
%% algorithm 3 with different q [estimated bound from paper sai19 theorem 1]

figure(1);
diffe  = zeros(1,3);
for q = 0:2
    [Uq,~,~,G] = AERandSVD(A,max,p,q);
    [stU,~] = angle_bounds(V,G,s,max,q);       
    sinthetau = subspace_angles(U(:,1:max),Uq); 
    semilogy(1:max,real(sinthetau),strcat('-'),...
            1:max,stU + 1.e-16,strcat('--'));
    diffe(q+1) = norm(stU + 1.e-16 - real(sinthetau));
    hold on
end
    
legend({'Computed q = 0', 'Estimate q = 0', ...
    'Computed q = 1', 'Estimate q = 1', ...
    'Computed q = 2', 'Estimate q = 2'},'Location','southeast');
ylabel('$\sin\theta_j$','Interpreter','LaTeX')
xlabel('Index')
savefig(figure(1),'Exp3alg3diffqthm1.fig')
diffe
%% algorithm 4 with different q [estimated bound from paper sai19 theorem 1]
%(I still use the same estimated bound since alg 4 is basically the same as alg 3 but eliminate the round off error when q is big)

figure(2);
diffe  = zeros(1,3);
for q = 0:2
    [Uq,~,~,G] = AEORandSVD(A,max,p,q);
    [stU,~] = angle_bounds(V,G,s,max,q);       
    sinthetau = subspace_angles(U(:,1:max),Uq); 
    semilogy(1:max,real(sinthetau),strcat('-'),...
            1:max,stU + 1.e-16,strcat('--'));
    diffe(q+1) = norm(stU + 1.e-16 - real(sinthetau));
    hold on
end
    
legend({'Computed q = 0', 'Estimate q = 0', ...
    'Computed q = 1', 'Estimate q = 1', ...
    'Computed q = 2', 'Estimate q = 2'},'Location','southeast');
ylabel('$\sin\theta_j$','Interpreter','LaTeX')
xlabel('Index')
savefig(figure(2),'Exp3alg4diffqthm1.fig')
diffe
%% algorithm 2,3,4 with q = 0 
% if q is 0 then basically 2,3,4 are the same thing!!
q = 0;
figure(3);
[Uq2,~,~,G2] = BasicRandSVD(A,max,p);      
sinthetau2 = subspace_angles(U(:,1:max),Uq2); 
semilogy(1:max,real(sinthetau2),'-');
hold on;
[Uq3,~,~,G3] = AERandSVD(A,max,p,q);     
sinthetau3 = subspace_angles(U(:,1:max),Uq3); 
semilogy(1:max,real(sinthetau3),'--');
[Uq4,~,~,G4] = AEORandSVD(A,max,p,q);       
sinthetau4 = subspace_angles(U(:,1:max),Uq4); 
semilogy(1:max,real(sinthetau4),'-.');

legend({'Computed Algorithm 2', ...
    'Computed Algorithm 3', ...
    'Computed Algorithm 4'},'Location','southeast');
ylabel('$\sin\theta_j$','Interpreter','LaTeX')
xlabel('Index')
savefig(figure(3),'Exp3alg234.fig')