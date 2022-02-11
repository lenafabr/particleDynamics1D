% PLOTSIMDENSITIES processes and plots densities of particles moving along
% tubular 1D domains
% requires *.snap.out file generated using partdynamics1D.exe
%% add required scripts to path
addpath('./tools/')
%% simulation output files
fnamestr = './examples/example.snap.out';
options = struct('readlastsnap',1);
%%
tic
			
[grouplist,tvals,domlen,ntrials] = readsnapshot(fnamestr,options);

% read particle positions
phagopos = [];
lysopos = [];
phagostates = [];
phtdeg = [];
phxfuse = [];
phacid = [];
phprot = [];
phnfuse = [];
lysointip = [];

for tc = 1:length(grouplist)
	pos = vertcat(grouplist(tc).particles(:).pos);
	states = vertcat(grouplist(tc).particles(:).state);
	ptypes = vertcat(grouplist(tc).particles(:).type);
	acid = vertcat(grouplist(tc).particles(:).acidity);
	prot = vertcat(grouplist(tc).particles(:).prot);
	prot = prot(:,2);
	nfuse = vertcat(grouplist(tc).particles(:).nfuse);
	tdeg = vertcat(grouplist(tc).particles(:).tdeg);
	xfuse = vertcat(grouplist(tc).particles(:).xfuse);
  intip = vertcat(grouplist(tc).particles(:).intip);
	
	phagopos = cat(1,phagopos,pos(ptypes==2));
	lysopos = cat(1,lysopos,pos(ptypes==1 & intip==0));
	phagostates = cat(1,phagostates,states(ptypes==2));
	phacid = cat(1,phacid,acid(ptypes==2));
	phprot = cat(1,phprot,prot(ptypes==2));
	phnfuse = cat(1,phnfuse,nfuse(ptypes==2));
	phtdeg = cat(1,phtdeg,tdeg(ptypes==2));
	phxfuse = cat(1,phxfuse,xfuse(ptypes==2));
  lysointip = cat(1,lysointip,intip(ptypes==1));
end
save('example_processed.mat');
toc
%% set parameters to discretize domain
vp = 1;
vy = 2.67;
Lreal = 1055;
dt = 1d-4; % simulation timestep
lyfrate = 100.001; % framerate for lysosomes (binwidth)
phfrate = (vy/vp)*lyfrate; % framerate for phagosomes (binwidth)
xp = 0:(phfrate*vp*dt):1; % possible positions for phagosomes
xp(end) = 1;
xl = 0:(lyfrate*vy*dt):1; % possible positions for lysosomes
xl(end) = 1;
%% plot organelle densities
clf
% organelle densities from sims
[nphago,~,bininds] = histcounts(phagopos,xp);
nfusehist = accumarray(bininds,phnfuse)'./nphago;
xpvals = 0.5*(xp(1:end-1)+xp(2:end));
dxp = (xpvals(2)-xpvals(1));

cphago = [215,38,156]/255;
clyso = [126,212,238]/255;
plot(xpvals*Lreal,fliplr(nphago)./ntrials./dxp./Lreal,'o-','Color',cphago,'MarkerFaceColor',cphago)

nfusedphago = histcounts(phagopos(phnfuse>0),xl);
nlyso = histcounts(lysopos,xl);
xlvals = 0.5*(xl(1:end-1)+xl(2:end));
dxl = (xlvals(2)-xlvals(1));

hold on
plot(xlvals*Lreal,fliplr(nlyso+nfusedphago)./ntrials./dxl./Lreal,'o-','Color',clyso,'MarkerFaceColor',clyso)
hold off

legend('AV density', 'lysosome density','FontSize',28)

plot_cleanup('Interpreter','tex','FontName','Arial','FontSize',28)
xlim([0,1050]);
ylabel('organelle density (\mum^{-1})')
xlabel('distance from distal end (\mum)')
%% plot fraction fused over space
plot(xpvals*Lreal,nfusehist,'bo-')
plot_cleanup('Interpreter','tex','FontName','Arial');
xlabel('position (\mum)')
ylabel('fraction of AVs fused')
ylim([0,1])
xlim([0,1060])

%% location of first fusion
binwidth = 100;
% xfusevals = (1-phxfuse(phxfuse>0 & phxfuse<1))*params.Lreal; % fused in the domain (excludes tip fusions)
xfusevals = (1-phxfuse(phxfuse>0))*Lreal; % fused anywhere
h = histogram(xfusevals,'BinWidth',binwidth,'Normalization','probability');
xlabel('distance from distal tip at the time of fusion (\mum)')
ylabel('fraction of fused AVs')
plot_cleanup('Interpreter','tex','FontName','Arial','FontSize',28)