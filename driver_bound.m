clear all

%% Plot controlled gap
m = 3000; 
n = 300;  
max = 35;
r = 15;
gap = [1,2,10];
figure(1);
lstyle = ["-","--",":"];

for i = 1:3
    A = controlledgap(m,n,r,gap(i));   
    s = svds(A,max+1);%subspace svd with first 35 singular value 
    semilogy(1:max,s(1:max),lstyle(i));
    hold on;
end
legend('Gap = 1', 'Gap = 2', 'Gap = 10');
xlabel('Index');
ylabel('Singular Values');
savefig(figure(1),'B1.fig')
%% general A for following test
m = 3000; 
n = 300; 
max = 35;
r = 15;
gap = 10;
A = controlledgap(m,n,r,gap); 
s = svds(A,max);
%% fixing p, varying k thm 2
p = 5;
k = [5,10,15,20,25,30];
L = zeros(1,6);
R = zeros(1,6);
for j = 1:6
    Q = FixedRank(A,k(j),p);
    L(j) = norm(A-Q*Q'*A);
    R(j) = (1+(4*sqrt(k(j)+p))/(p-1)*sqrt(min(m,n)))*s(k(j)+1);
end

figure(2);
h1 = semilogy(k,L,'r-');
hold on;
h2 = semilogy(k,R,'b--');
legend([h1(1), h2(1)],'Computed error', 'Estimated error');
xlabel('Values of k');
ylabel('Errors');
savefig(figure(2),'Bfixedp1.fig')
%% fixing k, varying p thm 2
k = 20;
p = [5,10,15,20,25];
L = zeros(1,5);
R = zeros(1,5);
for q = 1:5
    Q = FixedRank(A,k,p(q));
    L(q) = norm(A-Q*Q'*A);
    R(q) = (1+(4*sqrt(k+p(q)))/(p(q)-1)*sqrt(min(m,n)))*s(k+1);
end
figure(3);
h1 = semilogy(p,L,'r-');
hold on;
h2 = semilogy(p,R,'b--');

legend([h1(1), h2(1)],'Computed error', 'Estimated error');
xlabel('Values of p');
ylabel('Errors');
%axis mode
savefig(figure(3),'Bfixedk1.fig')
%% 100 runs with fixed p&k thm 2
k = 20;
p = 5;
bound = (1+(4*sqrt(k+p))/(p-1)*sqrt(min(m,n)))*s(k+1);
E = zeros(1,100);

for w = 1:100
     Q = FixedRank(A,k,p);
     E(w) = norm(A-Q*Q'*A);
end
averageE = sum(E)/100;

figure(4);
h1 = semilogy(1:w,E,'r-');
hold on;
h2 = yline(bound,'g--');
hold on;
h3 = yline(averageE,'b:');

legend([h1(1), h2(1), h3(1)],'Computed error for each run', 'Fixed Estimated error', 'Average Computed error over 100 runs');
xlabel('Runs')
ylabel('Errors')
axis([0 100 0 bound+10])

savefig(figure(4),'B100G1.fig')

%% fixing p,q, varying k thm 3
p = 5;
q = 1;
k = [5,10,15,20,25,30];
L = zeros(1,6);
R = zeros(1,6);
for j = 1:6
    Q = fixedrank_poweriteration(A,k(j),p,q);
    L(j) = norm(A-Q*Q'*A);
    R(j) = (1+ (4*sqrt((k(j)+p)*min(m,n)))/(p-1))^(1/(2*q+1))*s(k(j)+1);
end

figure(5);
h1 = semilogy(k,L,'r-');
hold on;
h2 = semilogy(k,R,'b--');
legend([h1(1), h2(1)],'Computed error', 'Estimated error');
xlabel('Values of k');
ylabel('Errors');
savefig(figure(5),'Bfixedpq2.fig')

%% fixing k,q, varying p thm 3
k = 20;
q = 1;
p = [5,10,15,20,25];
L = zeros(1,5);
R = zeros(1,5);
for j = 1:5
    Q = fixedrank_poweriteration(A,k,p(j),q);
    L(j) = norm(A-Q*Q'*A);
    R(j) = (1+ (4*sqrt((k+p(j))*min(m,n)))/(p(j)-1))^(1/(2*q+1))*s(k+1);
end
figure(6);
h1 = semilogy(p,L,'r-');
hold on;
h2 = semilogy(p,R,'b--');

legend([h1(1), h2(1)],'Computed error', 'Estimated error');
xlabel('Values of p');
ylabel('Errors');
%axis mode
savefig(figure(6),'Bfixedkq2.fig')

%% fixing k,p, varying q thm 3
k = 20;
p = 5;
q = [1,2];
L = zeros(1,2);
R = zeros(1,2);
for j = 1:2
    Q = fixedrank_poweriteration(A,k,p,q(j));
    L(j) = norm(A-Q*Q'*A);
    R(j) = (1+ (4*sqrt((k+p)*min(m,n)))/(p-1))^(1/(2*q(j)+1))*s(k+1);
end
figure(7);
h1 = semilogy(q,L,'r-');
hold on;
h2 = semilogy(q,R,'b--');

legend([h1(1), h2(1)],'Computed error', 'Estimated error');
xlabel('Values of q');
ylabel('Errors');
%axis mode
savefig(figure(7),'Bfixedkp2.fig')

%% 100 runs with fixed p,q&k thm 3
k = 20;
p = 5;
q = 1;
bound = (1+ (4*sqrt((k+p)*min(m,n)))/(p-1))^(1/(2*q+1))*s(k+1);
E = zeros(1,100);

for w = 1:100
     Q = fixedrank_poweriteration(A,k,p,q);
     E(w) = norm(A-Q*Q'*A);
end
averageE = sum(E)/100;

figure(8);
h1 = semilogy(1:w,E,'r-');
hold on;
h2 = yline(bound,'g--');
hold on;
h3 = yline(averageE,'b:');

legend([h1(1), h2(1), h3(1)],'Computed error for each run', 'Fixed Estimated error', 'Average Computed error over 100 runs');
xlabel('Runs')
ylabel('Errors')
axis([0 100 0 bound+1])

savefig(figure(8),'B100G2.fig')

%% fixing p,q, varying k thm 4
p = 5;
q = 1;
k = [5,10,15,20,25,30];
L = zeros(1,6);
R = zeros(1,6);
for j = 1:6
    Q = fixedrank_poweriterationO(A,k(j),p,q);
    L(j) = norm(A-Q*Q'*A);
    R(j) = ((1+sqrt(k(j)/(p-1)))+((exp(1)*sqrt(k(j)+p))/p)*sqrt(min(m,n)-k(j)))^(1/(2*q+1))*s(k(j)+1);
end

figure(9);
h1 = semilogy(k,L,'r-');
hold on;
h2 = semilogy(k,R,'b--');
legend([h1(1), h2(1)],'Computed error', 'Estimated error');
xlabel('Values of k');
ylabel('Errors');
savefig(figure(9),'Bfixedpq3.fig')

%% fixing k,q, varying p thm 4
k = 20;
q = 1;
p = [5,10,15,20,25];
L = zeros(1,5);
R = zeros(1,5);
for j = 1:5
    Q = fixedrank_poweriterationO(A,k,p(j),q);
    L(j) = norm(A-Q*Q'*A);
    R(j) = ((1+sqrt(k/(p(j)-1)))+((exp(1)*sqrt(k+p(j)))/p(j))*sqrt(min(m,n)-k))^(1/(2*q+1))*s(k+1);
end
figure(10);
h1 = semilogy(p,L,'r-');
hold on;
h2 = semilogy(p,R,'b--');

legend([h1(1), h2(1)],'Computed error', 'Estimated error');
xlabel('Values of p');
ylabel('Errors');
%axis mode
savefig(figure(10),'Bfixedkq3.fig')

%% fixing k,p, varying q thm 4
k = 20;
p = 5;
q = [1,2];
L = zeros(1,2);
R = zeros(1,2);
for j = 1:2
    Q = fixedrank_poweriterationO(A,k,p,q(j));
    L(j) = norm(A-Q*Q'*A);
    R(j) = ((1+sqrt(k/(p-1)))+((exp(1)*sqrt(k+p))/p)*sqrt(min(m,n)-k))^(1/(2*q(j)+1))*s(k+1);
end
figure(11);
h1 = semilogy(q,L,'r-');
hold on;
h2 = semilogy(q,R,'b--');

legend([h1(1), h2(1)],'Computed error', 'Estimated error');
xlabel('Values of q');
ylabel('Errors');
%axis mode
savefig(figure(11),'Bfixedkp3.fig')

%% 100 runs with fixed p,q&k thm 4
k = 20;
p = 5;
q = 1;
bound = ((1+sqrt(k/(p-1)))+((exp(1)*sqrt(k+p))/p)*sqrt(min(m,n)-k))^(1/(2*q+1))*s(k+1);
E = zeros(1,100);

for w = 1:100
     Q = fixedrank_poweriterationO(A,k,p,q);
     E(w) = norm(A-Q*Q'*A);
end
averageE = sum(E)/100;

figure(12);
h1 = semilogy(1:w,E,'r-');
hold on;
h2 = yline(bound,'g--');
hold on;
h3 = yline(averageE,'b:');

legend([h1(1), h2(1), h3(1)],'Computed error for each run', 'Fixed Estimated error', 'Average Computed error over 100 runs');
xlabel('Runs')
ylabel('Errors')
axis([0 100 0 bound+1])

savefig(figure(12),'B100G3.fig')