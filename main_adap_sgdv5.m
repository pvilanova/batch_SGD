close all;
clearvars;

        x = [20, 50]/10;
        d = numel(x);
        
        ka = 100;
        H = @(t) eye(d)*(1-t) + [2*ka 1/2; 1/2 1]*t;
        f = @(x,t) x*H(t)*x'/2 - ones(1,d)*x';

        x_so = (H(1/2)\ones(d,1))';
        L = max(eig(H(1/2)));
        mu = min(eig(H(1/2)));

%         nu = @(t) 2/(1+L*t);

        md = 1/2;
        get_samp = @(M) unifrnd(0,1,M,1); 
        
        x_s = sym('x', [1 d], 'Real');
        t_s = sym('t', 'Real');
        tmp = jacobian(f(x_s,t_s),x_s);
        gf = matlabFunction(tmp, 'Vars', {[x_s], t_s});

       cov_matrix = H(1).^2;
tr_sig = trace(cov_matrix);

% v_ep = [0.1 0.5 1 5.908942376];
% v_nu = [0.08660254 0.433012702 0.866025404 5.84];
% v_vartheta = [0.05 0.25 0.5 0.9];

v_ep = [0.5 0.9];
v_vartheta = v_ep/2;
v_nu = sqrt(v_ep.^2-v_vartheta.^2);

v_par = [v_ep; v_nu; v_vartheta];

%v_its = [40 20 30 400];
v_its = [20 20];
    
%j = 1;

%marks = {'+','o','*','.','x','s','d','^'};
line = {'-',':','--','-.'};

for j=1:size(v_par,2)

ep = v_par(1,j);
nu = v_par(2,j);
vartheta = v_par(3,j);

%for ep = v_ep
% for nu = v_nu
% for vartheta = v_vartheta

%vartheta = ep^2 - nu^2;
%nu = sqrt(ep^2 - vartheta^2);
%nu = 5.84;

its = v_its(j);

reps = 1E4;
gapn = zeros(reps,its);
gapd = zeros(reps,its);
errn = zeros(reps,its);     % verification purposes
errd = zeros(reps,its);     %
v_bn = zeros(reps,its);
v_bd = zeros(reps,its);

tic

for kk=1:reps
    
rng(kk)

xn = x;
for k=1:its
    grad_norm = norm(gf(xn(end,:),md));
    
    b = ceil(tr_sig/(ep^2*grad_norm^2));
    
    thetas = get_samp(b);
    
    v_bn(kk,k) = b;
    step = 2/(mu+L)/(1+ep^2);
    xn = [xn; xn(end,:) - step*mean(gf(xn(end,:), thetas),1)];
end

rng(kk)

xd = x;
for k=1:its
    gr = gf(xd(end,:),md);
    grad_norm = norm(gr);

    e_grad = gr/grad_norm;
    
    sig_grad = e_grad*cov_matrix*e_grad';
    sig_perp = tr_sig - sig_grad;
    b_grad = ceil(sig_grad/(vartheta^2*grad_norm^2));
    b_perp = ceil(sig_perp/(nu^2*grad_norm^2));
    b = min(b_grad,b_perp);
    
    thetas = get_samp(b);
    errd(kk,k) = mean(((gf(xd(end,:), thetas)-gr)*e_grad').^2);
    
    v_bd(kk,k) = b;
    step = 2/(mu+L)/(1+nu^2+vartheta^2);
    xd = [xd; xd(end,:) - step*mean(gf(xd(end,:), thetas),1)];
end

for k=1:its
    gapn(kk,k) = abs(f(xn(k,:),md) - f(x_so,md));
    gapd(kk,k) = abs(f(xd(k,:),md) - f(x_so,md));
end

end

toc

qgapn = zeros(2,its);
qgapd = zeros(2,its);
for k=1:its
    tmpn = bootstrp(reps, @mean, gapn(:,k));
    tmpd = bootstrp(reps, @mean, gapd(:,k));
    qgapn(:,k) = [quantile(tmpn,0.025)  quantile(tmpn,0.975)];
    qgapd(:,k) = [quantile(tmpd,0.025)  quantile(tmpd,0.975)];
end

qv_bn = zeros(2,its);
qv_bd = zeros(2,its);
for k=1:its
    tmpn = bootstrp(reps, @mean, v_bn(:,k));
    tmpd = bootstrp(reps, @mean, v_bd(:,k));
    qv_bn(:,k) = [quantile(tmpn,0.025)  quantile(tmpn,0.975)];
    qv_bd(:,k) = [quantile(tmpd,0.025)  quantile(tmpd,0.975)];
end

% tit_str = sprintf('epsilon = %g, theta = %g, nu = %g',ep,vartheta,nu);

titd = sprintf('\\epsilon=%g, \\vartheta=%g, \\nu=%g',ep,vartheta,nu);
titn = sprintf('\\epsilon=%g',ep);

rate = ((L-mu)^2/(L+mu)^2+ep^2)/(1+ep^2);
bound = (L/2)*rate.^[1:its]*norm(x)^2;

lin = line{1 + mod(j-1, numel(line))};
figure(1);
hold on;
plot(bound,['k' lin],'DisplayName',[titn ' rate']);
plot(mean(gapn,1),['r' lin],'DisplayName',[titn ' norm']);
ciplot(mean(qgapn(1,:),1),mean(qgapn(2,:),1),1:its,'r');
% semilogy(mean(qgapn(1,:),1),['r:' marks{mod(j,numel(marks))}],'DisplayName','95% low CI norm');
% semilogy(mean(qgapn(2,:),1),['r:' marks{mod(j,numel(marks))}],'DisplayName','95% high CI norm');
plot(mean(gapd,1),['b' lin],'DisplayName',[titd ' inner/orth']);
ciplot(mean(qgapd(1,:),1),mean(qgapd(2,:),1),1:its,'b');
% semilogy(mean(qgapd(1,:),1),['b:' marks{mod(j,numel(marks))}],'DisplayName','95% low CI inner/orth');
% semilogy(mean(qgapd(2,:),1),['b:' marks{mod(j,numel(marks))}],'DisplayName','95% high CI inner/orth');
hold off;
xlabel('iteration');
ylabel('average gap');

figure(2);
hold on;
plot(mean(v_bn,1),['r' lin],'DisplayName',[titn ' norm']);
ciplot(mean(qv_bn(1,:),1),mean(qv_bn(2,:),1),1:its,'r');
%semilogy(mean(qv_bn(1,:),1),['r:' marks{mod(j,numel(marks))}],'DisplayName','95% low CI norm');
%semilogy(mean(qv_bn(2,:),1),['r:' marks{mod(j,numel(marks))}],'DisplayName','95% high CI norm');
plot(mean(v_bd,1),['b' lin],'DisplayName',[titd ' inner/orth']);
ciplot(mean(qv_bd(1,:),1),mean(qv_bd(2,:),1),1:its,'b');
%semilogy(mean(qv_bd(1,:),1),['b:' marks{mod(j,numel(marks))}],'DisplayName','95% low CI norm');
%semilogy(mean(qv_bd(2,:),1),['b:' marks{mod(j,numel(marks))}],'DisplayName','95% high CI norm');
hold off;
xlabel('iteration');
ylabel('average b per level');

figure(3);
hold on;
ciplot(mean(qgapn(1,:),1),mean(qgapn(2,:),1),cumsum(mean(v_bn,1)),'r');
plot(cumsum(mean(v_bn,1)),mean(gapn,1),['r' lin],'DisplayName',[titn ' norm']);
% loglog(cumsum(mean(qv_bn(1,:),1)),mean(qgapn(1,:),1),['r:' marks{mod(j,numel(marks))}],'DisplayName','95% low CI norm');
% loglog(cumsum(mean(qv_bn(2,:),1)),mean(qgapn(2,:),1),['r:' marks{mod(j,numel(marks))}],'DisplayName','95% high CI norm');
ciplot(mean(qgapd(1,:),1),mean(qgapd(2,:),1),cumsum(mean(v_bd,1)),'b');
plot(cumsum(mean(v_bd,1)),mean(gapd,1),['b' lin],'DisplayName',[titd, ' inner/orth']);
% loglog(cumsum(mean(qv_bd(1,:),1)),mean(qgapd(1,:),1),['b:' marks{mod(j,numel(marks))}],'DisplayName','95% low CI norm');
% loglog(cumsum(mean(qv_bd(2,:),1)),mean(qgapd(2,:),1),['b:' marks{mod(j,numel(marks))}],'DisplayName','95% high CI norm');
hold off;
xlabel('average total b','FontSize',18);
ylabel('average gap','FontSize',18);
title('Optimal gap vs. total samples','FontSize',20);

bound = zeros(1,its);
for k=1:its
    bound(k) = norm(gf(xn(k,:),md))^2;
end
figure(4);
hold on;
plot(ep^2*bound,['k' lin],'DisplayName',[titn ' bound']);
plot(mean(errn,1),['r' lin],'DisplayName',[titn ' norm']);
plot(mean(errd,1),['b' lin],'DisplayName',[titd ' inner/orth']);
hold off;
xlabel('iteration');
ylabel('norm test per level');
%j = j+1;

end
% end
% end

figure(1);
legend('location','best');
set(gca,'yscale','log');
savefig(gcf,'its_gap');
print(gcf,'its_gap.eps','-depsc');

figure(2);
legend('location','best');
set(gca,'yscale','log');
savefig(gcf,'its_b');
print(gcf,'its_b.eps','-depsc');

figure(3);
legend('location','best');
set(gca,'yscale','log');
set(gca,'xscale','log');
savefig(gcf,'b_gap');
print(gcf,'b_gap.eps','-depsc');

figure(4);
legend('location','best');
set(gca,'yscale','log');
savefig(gcf,'its_err');
print(gcf,'its_err.eps','-depsc');

