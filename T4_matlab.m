clc; close all; clear
warning off
spiralPitch = 1.7; 
kValue = spiralPitch/2/pi; 
mainLength = 341e-2;
holeDist1 = mainLength - 27.5e-2*2; 
secLength = 220e-2;
holeDist2 = secLength - 27.5e-2*2; 
initVelocity = 1; 

angleSet = 5*2*pi:-0.01:0*pi;
radius = kValue*angleSet;
xCoord = radius.*cos(angleSet);
yCoord = radius.*sin(angleSet);
figure(1)
set(gcf, 'Position', [200 200 600 600]);
colorSet1 = rand(1,3);
plot(xCoord, yCoord, '-', 'Color', colorSet1, 'LineWidth', 1.3)
axis equal
grid on
xlabel('x')
ylabel('y')
set(gca, 'FontSize', 18) 
hold on

angleSet = angleSet - pi;
radius2 = kValue*(angleSet + pi); 
xCoord2 = radius2.*cos(angleSet);
yCoord2 = radius2.*sin(angleSet);
colorSet2 = rand(1,3);
plot(xCoord2, yCoord2, 'm', 'Color', colorSet2, 'LineWidth', 1.3)
turnRadius = 4.5; 
xTurn = turnRadius*cos(angleSet);
yTurn = turnRadius*sin(angleSet);
colorSet3 = sort(rand(1,3));
plot(xTurn, yTurn, 'Color', colorSet3, 'LineWidth', 2)

entryTheta = turnRadius/kValue;
exitTheta = turnRadius/kValue - pi;   
slopeValue = (kValue*sin(entryTheta) + turnRadius*cos(entryTheta)) / (kValue*cos(entryTheta) - turnRadius*sin(entryTheta)); 
maxTheta1 = atan(-1/slopeValue) + pi; 
angleEquality = atan(tan(entryTheta)) + pi - maxTheta1; 
totalRadius = turnRadius / cos(angleEquality); 
smallRadius = totalRadius / 3; 
largeRadius = smallRadius * 2; 
phiValue = 2 * angleEquality; 
arcLength1 = largeRadius * (pi - phiValue); 
arcLength2 = smallRadius * (pi - phiValue); 
minTheta1 = maxTheta1 - arcLength1 / largeRadius;  
minTheta2 = minTheta1 - pi; 
maxTheta2 = minTheta2 + arcLength2 / smallRadius; 
centerX1 = turnRadius*cos(entryTheta) + largeRadius*cos(maxTheta1 - pi);
centerY1 = turnRadius*sin(entryTheta) + largeRadius*sin(maxTheta1 - pi); 

centerX2 = turnRadius*cos(exitTheta) - smallRadius*cos(maxTheta2);
centerY2 = turnRadius*sin(exitTheta) - smallRadius*sin(maxTheta2); 
figure(1)
hold on
plot(centerX1 + largeRadius*cos(linspace(minTheta1, maxTheta1, 50)), centerY1 + largeRadius*sin(linspace(minTheta1, maxTheta1, 50)), 'r', 'LineWidth', 2)
plot(centerX1, centerY1, 'r*')
plot(centerX2 + smallRadius*cos(linspace(minTheta2, maxTheta2, 50)), centerY2 + smallRadius*sin(linspace(minTheta2, maxTheta2, 50)), 'b', 'LineWidth', 2)
plot(centerX2, centerY2, 'b*')
axis equal

diffTheta = @(t,theta) 1./(kValue*sqrt(1+theta.^2));
startTheta = entryTheta; 
timeStep = 0.1; 
timeStep = 1/randi([5 20]); 
timeSpan = 0:timeStep:100; 
[timeVals, thetaVals] = ode45(diffTheta, timeSpan, startTheta); 
xFirst = kValue*thetaVals.*cos(thetaVals);
yFirst = kValue*thetaVals.*sin(thetaVals);
timeReversed = timeVals(end:-1:1);
xMatrix = zeros(224, 200/timeStep+1); 
yMatrix = zeros(224, 200/timeStep+1); 
thetaMatrix = zeros(224, 200/timeStep+1);
xMatrix(1, 1:length(xFirst)) = xFirst(end:-1:1); 
yMatrix(1, 1:length(yFirst)) = yFirst(end:-1:1); 
thetaMatrix(1, 1:length(thetaVals)) = thetaVals(end:-1:1); 
timeArc1 = timeStep:timeStep:arcLength1; 
thetaArc1 = -timeArc1 / largeRadius + maxTheta1;
thetaMatrix(1, length(thetaVals)+(1:length(timeArc1))) = thetaArc1;
xMatrix(1, length(xFirst)+(1:length(timeArc1))) = largeRadius*cos(thetaArc1) + centerX1; 
yMatrix(1, length(yFirst)+(1:length(timeArc1))) = largeRadius*sin(thetaArc1) + centerY1; 
timeArc2 = timeArc1(end)+timeStep:timeStep:arcLength1+arcLength2; 
thetaArc2 = (timeArc2 - arcLength1) / smallRadius + minTheta1 - pi; 
thetaMatrix(1, length(thetaVals) + length(timeArc1) + (1:length(timeArc2))) = thetaArc2;
xMatrix(1, length(xFirst) + length(timeArc1) + (1:length(timeArc2))) = smallRadius*cos(thetaArc2) + centerX2; 
yMatrix(1, length(yFirst) + length(timeArc1) + (1:length(timeArc2))) = smallRadius*sin(thetaArc2) + centerY2; 
diffTheta = @(t,theta) 1./(kValue*sqrt(1+(theta+pi).^2)); 
startTheta = exitTheta; 
timeSpan = timeArc2(end)+timeStep:timeStep:100; 
[timeVals, theta2Vals] = ode45(diffTheta, timeSpan, startTheta); 
xSecond = kValue*(theta2Vals+pi).*cos(theta2Vals);
ySecond = kValue*(theta2Vals+pi).*sin(theta2Vals);
thetaMatrix(1, length(thetaVals) + length(timeArc1) + length(timeArc2) + 1:end) = theta2Vals;
xMatrix(1, length(xFirst) + length(timeArc1) + length(timeArc2) + 1:end) = xSecond;
yMatrix(1, length(yFirst) + length(timeArc1) + length(timeArc2) + 1:end) = ySecond; 
figure(3)
set(gcf, 'Position', [300 300 600 600])

for i=1:3:length(thetaMatrix(1,:))
    title({['t=', num2str((i-1)*timeStep)], 'VX公众号Matlab techniques出品', '头部把手中心的轨迹(-100到100s)'})
    plot(xMatrix(1,i), yMatrix(1,i), 'Marker', 'o', 'MarkerSize', 3, 'MarkerFaceColor', 'r')
    hold on
    axis equal
    axis([-10 10 -10 10])
    grid on
    drawnow
end

progressBar = waitbar(0, '计算开始...');
totalNumbers = 223; 
totalTimeSpan = -100:timeStep:100;
for jj = (1:length(totalTimeSpan))
    currentIndex = round((totalTimeSpan(1)+100)/timeStep)+jj;
    if totalTimeSpan(jj) <= 0 
        for i = 2:totalNumbers+1 
            dist = holeDist1*(i <= 2) + holeDist2*(i > 2); 
            nextTheta = solve_theta1(spiralPitch, xMatrix(i-1, currentIndex), yMatrix(i-1, currentIndex), thetaMatrix(i-1, currentIndex), dist); 
            thetaMatrix(i, currentIndex) = nextTheta;
            xMatrix(i, currentIndex) = kValue*nextTheta*cos(nextTheta);
            yMatrix(i, currentIndex) = kValue*nextTheta*sin(nextTheta);
        end
    elseif totalTimeSpan(jj) > 0 && totalTimeSpan(jj) <= arcLength1 
        currentFlag = 2;
        for i = 2:totalNumbers+1
            dist = holeDist1*(i <= 2) + holeDist2*(i > 2); 
            if currentFlag == 2 
                [nextX, nextY, nextTheta, currentFlag] = solve_point_2_1(spiralPitch, xMatrix(i-1, currentIndex), yMatrix(i-1, currentIndex), thetaMatrix(i-1, currentIndex), dist, largeRadius, centerX1, centerY1, maxTheta1);
                thetaMatrix(i, currentIndex) = nextTheta;
                xMatrix(i, currentIndex) = nextX;
                yMatrix(i, currentIndex) = nextY;
            else  
                nextTheta = solve_theta1(spiralPitch, xMatrix(i-1, currentIndex), yMatrix(i-1, currentIndex), thetaMatrix(i-1, currentIndex), dist); 
                thetaMatrix(i, currentIndex) = nextTheta;
                xMatrix(i, currentIndex) = kValue*nextTheta*cos(nextTheta);
                yMatrix(i, currentIndex) = kValue*nextTheta*sin(nextTheta);
            end
        end
    elseif totalTimeSpan(jj) > arcLength1 && totalTimeSpan(jj) <= arcLength1 + arcLength2 
        currentFlag = 3;
        for i = 2:totalNumbers+1
            dist = holeDist1*(i <= 2) + holeDist2*(i > 2); 
            if currentFlag == 3 
                [nextX, nextY, nextTheta, currentFlag] = solve_point_3_2(xMatrix(i-1, currentIndex), yMatrix(i-1, currentIndex), thetaMatrix(i-1, currentIndex), dist, largeRadius, centerX1, centerY1, smallRadius, centerX2, centerY2, minTheta2);
                thetaMatrix(i, currentIndex) = nextTheta;
                xMatrix(i, currentIndex) = nextX;
                yMatrix(i, currentIndex) = nextY;
            elseif currentFlag == 2 
                [nextX, nextY, nextTheta, currentFlag] = solve_point_2_1(spiralPitch, xMatrix(i-1, currentIndex), yMatrix(i-1, currentIndex), thetaMatrix(i-1, currentIndex), dist, largeRadius, centerX1, centerY1, maxTheta1);
                thetaMatrix(i, currentIndex) = nextTheta;
                xMatrix(i, currentIndex) = nextX;
                yMatrix(i, currentIndex) = nextY;
            else  
                nextTheta = solve_theta1(spiralPitch, xMatrix(i-1, currentIndex), yMatrix(i-1, currentIndex), thetaMatrix(i-1, currentIndex), dist); 
                thetaMatrix(i, currentIndex) = nextTheta;
                xMatrix(i, currentIndex) = kValue*nextTheta*cos(nextTheta);
                yMatrix(i, currentIndex) = kValue*nextTheta*sin(nextTheta);
            end
        end
    else  
        currentFlag = 4;
        for i = 2:totalNumbers+1
            dist = holeDist1*(i <= 2) + holeDist2*(i > 2); 
            if currentFlag == 4 
                [nextX, nextY, nextTheta, currentFlag] = solve_point_4_3(spiralPitch, xMatrix(i-1, currentIndex), yMatrix(i-1, currentIndex), thetaMatrix(i-1, currentIndex), dist, smallRadius, centerX2, centerY2, maxTheta2);
                thetaMatrix(i, currentIndex) = nextTheta;
                xMatrix(i, currentIndex) = nextX;
                yMatrix(i, currentIndex) = nextY;
            elseif currentFlag == 3 
                [nextX, nextY, nextTheta, currentFlag] = solve_point_3_2(xMatrix(i-1, currentIndex), yMatrix(i-1, currentIndex), thetaMatrix(i-1, currentIndex), dist, largeRadius, centerX1, centerY1, smallRadius, centerX2, centerY2, minTheta2);
                thetaMatrix(i, currentIndex) = nextTheta;
                xMatrix(i, currentIndex) = nextX;
                yMatrix(i, currentIndex) = nextY;
            elseif currentFlag == 2 
                [nextX, nextY, nextTheta, currentFlag] = solve_point_2_1(spiralPitch, xMatrix(i-1, currentIndex), yMatrix(i-1, currentIndex), thetaMatrix(i-1, currentIndex), dist, largeRadius, centerX1, centerY1, maxTheta1);
                thetaMatrix(i, currentIndex) = nextTheta;
                xMatrix(i, currentIndex) = nextX;
                yMatrix(i, currentIndex) = nextY;
            else  
                nextTheta = solve_theta1(spiralPitch, xMatrix(i-1, currentIndex), yMatrix(i-1, currentIndex), thetaMatrix(i-1, currentIndex), dist); 
                thetaMatrix(i, currentIndex) = nextTheta;
                xMatrix(i, currentIndex) = kValue*nextTheta*cos(nextTheta);
                yMatrix(i, currentIndex) = kValue*nextTheta*sin(nextTheta);
            end
        end
    end
    waitbar(jj/length(totalTimeSpan), progressBar, '已经完成...')
end
close(progressBar)
figure(100)
clf;
set(gcf, 'Position', [200 200 600 600])
for j = (1:2:length(totalTimeSpan)) + round((totalTimeSpan(1) + 100) / timeStep)
    plot(xMatrix(:, j), yMatrix(:, j), 'k-', 'LineWidth', 1.2, 'Marker', 'o', 'MarkerSize', 6, 'MarkerFaceColor', 'r')
    title({['t=', num2str(timeStep * (j - 1) - 100)], 'VX公众号Matlab techniques出品', '盘入-掉头-盘出的轨迹'})
    axis equal
    grid on
    xlabel('x')
    ylabel('y')
    axis([-15 15 -15 15])
    drawnow
end

xVelocity = zeros(size(xMatrix)); 
yVelocity = xVelocity;
xVelocity(:, 1) = (xMatrix(:, 2) - xMatrix(:, 1)) / timeStep; 
xVelocity(:, end) = (xMatrix(:, end) - xMatrix(:, end - 1)) / timeStep; 
xVelocity(:, 2:end - 1) = (xMatrix(:, 3:end) - xMatrix(:, 1:end - 2)) / 2 / timeStep; 
yVelocity(:, 1) = (yMatrix(:, 2) - yMatrix(:, 1)) / timeStep; 
yVelocity(:, end) = (yMatrix(:, end) - yMatrix(:, end - 1)) / timeStep; 
yVelocity(:, 2:end - 1) = (yMatrix(:, 3:end) - yMatrix(:, 1:end - 2)) / 2 / timeStep; 
totalVelocity = sqrt(xVelocity.^2 + yVelocity.^2); 

figure
plot(totalTimeSpan, totalVelocity(1, :), 'b.', 'LineWidth', 1.3)
ylim([0 1.1])
xlabel('时间')
ylabel('头把手的速度')
title('验证数值计算得到的速度公众号：Matlab techniques')
format long
timeIndex = 1 / timeStep; 
timePoints = 1:timeIndex:length(totalTimeSpan); 
locationData = zeros(2 * (totalNumbers + 1), length(timePoints)); 
locationData(1:2:end, :) = round(xMatrix(:, timePoints), 6); 
locationData(2:2:end, :) = round(yMatrix(:, timePoints), 6); 
velocityData = round(totalVelocity(:, timePoints), 6); 
dataFilename = 'result44_test.xlsx';
dataSheet = 1;
excelRange = 'B2';
xlswrite(dataFilename, locationData, dataSheet, excelRange)  
dataSheet = 2;
excelRange = 'B2';
xlswrite(dataFilename, velocityData, dataSheet, excelRange)  

timeInterval = 50 / timeStep; 
timePoints2 = 1:timeInterval:length(totalTimeSpan); 
rowIndices = [1 2:50:totalNumbers 224];
locationData2 = zeros(2 * length(rowIndices), length(timePoints2)); 
locationData2(1:2:end, :) = round(xMatrix(rowIndices, timePoints2), 6); 
locationData2(2:2:end, :) = round(yMatrix(rowIndices, timePoints2), 6) 
velocityData2 = round(totalVelocity(rowIndices, timePoints2), 6)