data = load('data/fingerforcedata.txt');
% Take a look (first four columns are the force)
forces = data(:,1:4);
% Sample rate is 200 Hz
time = (1:size(forces,1)) .* 1/200;
plot(time,forces);
hold on;
sumforce = sum(forces,2);
plot(time,sumforce,'k');
xlabel('time (s)');
ylabel('Force (N)');
legend('Index','Middle','Ring','Little','Total');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[peaks,peakLocations] = ...
  findpeaks(sumforce,'minpeakheight',20,...
                     'minpeakdistance',100);
[troughs,troughLocations] = ...
  findpeaks(-sumforce,'minpeakheight',-10,...
                      'minpeakdistance',100);
figure;
plot(time,sumforce);
hold on;
plot(time(peakLocations),...
  sumforce(peakLocations),'r*');
plot(time(troughLocations),...
  sumforce(troughLocations),'g*');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=2:numel(troughLocations)-1
 movementsUp{k-1} = forces(troughLocations(k):peakLocations(k),:);
 timeUp{k-1} = time(troughLocations(k):peakLocations(k));
  for m=1:4
   movementsUpCombined(k-1,m,:) = ...
    interp1(timeUp{k-1},movementsUp{k-1}(:,m),...
    linspace(timeUp{k-1}(1),timeUp{k-1}(end),100));
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=2:numel(troughLocations)-1
 movementsDown{k-1} = forces(peakLocations(k):troughLocations(k+1),:);
 timeDown{k-1} = time(peakLocations(k):troughLocations(k+1));
 for m=1:4
  movementsDownCombined(k-1,m,:) = ...
   interp1(timeDown{k-1},movementsDown{k-1}(:,m),...
   linspace(timeDown{k-1}(1),timeDown{k-1}(end),100));
 end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms fi fm fr fl;
ftot = fi + fm + fr + fl;
J(1,1) = diff(ftot,fi);
J(1,2) = diff(ftot,fm);
J(1,3) = diff(ftot,fr);
J(1,4) = diff(ftot,fl);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
J = [1 1 1 1];
nullspace = null(J);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numrepetitions = size(movementsUpCombined,1);
meanMovementsUp = mean(movementsUpCombined);
meanFreeMovementsUp = movementsUpCombined - repmat(meanMovementsUp,numrepetitions,1);

meanMovementsDown = mean(movementsDownCombined);
meanFreeMovementsDown = movementsDownCombined - repmat(meanMovementsDown,numrepetitions,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t=1:100
  for k=1:size(nullspace,2)
    goodvariancesUp(t,k,:,:) = (nullspace(:,k)' * meanFreeMovementsUp(:,:,t)')' * nullspace(:,k)';
  end
  goodvariancesCombinedUp(:,:,t) = sum(goodvariancesUp(t,:,:,:));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vUCM_Up = squeeze(...
 sum(sum(goodvariancesCombinedUp.^2,2),1)) ...
 ./ ((4-1)*numrepetitions);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
badvariancesCombinedUp = meanFreeMovementsUp - goodvariancesCombinedUp;
vORT_Up = squeeze(...
  sum(sum(badvariancesCombinedUp.^2,2),1)) ...
  ./ (1*numrepetitions);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
deltaV_Up = (vUCM_Up-vORT_Up)./vORT_Up;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numrepetitions = size(movementsDownCombined,1);
meanMovementsDown = mean(movementsDownCombined);
meanFreeMovementsDown = movementsDownCombined - repmat(meanMovementsDown,numrepetitions,1);

meanMovementsDown = mean(movementsDownCombined);
meanFreeMovementsDown = movementsDownCombined - repmat(meanMovementsDown,numrepetitions,1);
for t=1:100
  for k=1:size(nullspace,2)
    goodvariancesDown(t,k,:,:) = (nullspace(:,k)' * meanFreeMovementsDown(:,:,t)')' * nullspace(:,k)';
  end
  goodvariancesCombinedDown(:,:,t) = sum(goodvariancesDown(t,:,:,:));
end
vUCM_Down = squeeze(sum(sum(goodvariancesCombinedDown.^2,2),1)) ./ ...
  ((4-1)*numrepetitions);
badvariancesCombinedDown = meanFreeMovementsDown - goodvariancesCombinedDown;
vORT_Down = squeeze(sum(sum(badvariancesCombinedDown.^2,2),1)) ...
  ./ (1*numrepetitions);
deltaV_Down = (vUCM_Down-vORT_Down)./vORT_Down;
%%
figure;
subplot(2,3,1);
plot(linspace(0,1,100),vUCM_Up);
ylabel('v_{UCM} Up');
xlabel('Normalized time');
subplot(2,3,2);
plot(linspace(0,1,100),vORT_Up);
ylabel('v_{ORT} Up');
xlabel('Normalized time');
subplot(2,3,3);
plot(linspace(0,1,100),deltaV_Up);
ylabel('\Delta V Up');
xlabel('Normalized time');

subplot(2,3,4);
plot(linspace(0,1,100),vUCM_Down);
ylabel('v_{UCM} Down');
xlabel('Normalized time');
subplot(2,3,5);
plot(linspace(0,1,100),vORT_Down);
ylabel('v_{ORT} Down');
xlabel('Normalized time');
subplot(2,3,6);
plot(linspace(0,1,100),deltaV_Down);
ylabel('\Delta V Down');
xlabel('Normalized time');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
