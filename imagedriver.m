%% read in matrix A and see it as a pic
%A = readmatrix('sunflower.txt');
[A, map] = imread('KSsunflower.png');
A = im2double(A);
[m,n] = size(A);
figure(1);
%image(A);

imshow(A);
title(['The actual figure with matrix A in size ',num2str(m),'x',num2str(n),' with rank ',num2str(rank(A))]);
saveas(figure(1),'Sunflower.png')

index = [10,50,100,400,800]
%% Test Alg 2
figure(2);
set(gcf, 'WindowState', 'maximized');
p = 10;
for i = 1:5
    k = index(i);
    tic;
    [U,D,V] = BasicRandSVD(A,k,p);
    t = toc;
    approxA = U*D*V';
    subplot(1,5,i)
    imshow(approxA);
    pbaspect([10 8 1])
    title(['p = ',num2str(p),', k = ',num2str(k),' in ',num2str(t),' time']);
end
%saveas(figure(2),'ALG2P10.png')

figure(3);
set(gcf, 'WindowState', 'maximized');
p = 0;
for i = 1:5
    k = index(i);
    tic
    [U,D,V] = BasicRandSVD(A,k,p);
    t = toc;
    approxA = U*D*V';
    subplot(1,5,i)
    imshow(approxA);
    pbaspect([10 8 1])
    title(['p = ',num2str(p),', k = ',num2str(k),' in ',num2str(t),' time']);
end
saveas(figure(3),'ALG2P0.png')

%% Test Alg 3
figure(4);
set(gcf, 'WindowState', 'maximized');
p = 10;
q = 1;
for i = 1:5
    k = index(i);
    tic;
    [U,D,V] = AERandSVD(A,k,p,q);
    t = toc;
    approxA = U*D*V';
    subplot(1,5,i)
    imshow(approxA);
    pbaspect([10 8 1])
    title(['q = ',num2str(q),', k = ',num2str(k),' in ',num2str(t),' time']);
end
saveas(figure(4),'ALG3Q1.png')

figure(5);
set(gcf, 'WindowState', 'maximized');
p = 10;
q = 2;
for i = 1:5
    k = index(i);
    tic
    [U,D,V] = BasicRandSVD(A,k,p);
    t = toc;
    approxA = U*D*V';
    subplot(1,5,i)
    imshow(approxA);
    pbaspect([10 8 1])
    title(['q = ',num2str(q),', k = ',num2str(k),' in ',num2str(t),' time']);
end
saveas(figure(5),'ALG3Q2.png')
%% Test Alg 4
figure(6);
set(gcf, 'WindowState', 'maximized');
p = 10;
q = 1;
for i = 1:5
    k = index(i);
    tic;
    [U,D,V] = AEORandSVD(A,k,p,q);
    t = toc;
    approxA = U*D*V';
    subplot(1,5,i)
    imshow(approxA);
    pbaspect([10 8 1])
    title(['q = ',num2str(q),', k = ',num2str(k),' in ',num2str(t),' time']);
end
saveas(figure(6),'ALG4Q1.png')

figure(7);
set(gcf, 'WindowState', 'maximized');
p = 10;
q = 2;
for i = 1:5
    k = index(i);
    tic
    [U,D,V] = BasicRandSVD(A,k,p);
    t = toc;
    approxA = U*D*V';
    subplot(1,5,i)
    imshow(approxA);
    pbaspect([10 8 1])
    title(['q = ',num2str(q),', k = ',num2str(k),' in ',num2str(t),' time']);
end
saveas(figure(7),'ALG4Q2.png')

%% Test Alg 6
figure(8);
set(gcf, 'WindowState', 'maximized');
p = 10;
for i = 1:5
    k = index(i);
    tic
    [U,D,V] = SPRandSVD(A,k,p);
    t = toc;
    approxA = U*D*V';
    subplot(1,5,i)
    imshow(approxA);
    pbaspect([10 8 1])
    title(['p = ',num2str(p),', k = ',num2str(k),' in ',num2str(t),' time']);
end
saveas(figure(8),'ALG6P10.png')

figure(9);
set(gcf, 'WindowState', 'maximized');
p = 0;
for i = 1:5
    k = index(i);
    tic
    [U,D,V] = SPRandSVD(A,k,p);
    t = toc;
    approxA = U*D*V';
    subplot(1,5,i)
    imshow(approxA);
    pbaspect([10 8 1])
    title(['p = ',num2str(p),', k = ',num2str(k),' in ',num2str(t),' time']);
end
saveas(figure(8),'ALG6P0.png')
%% 