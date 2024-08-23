function [shiftDf] = timePermute(df)


% Circular shift each cell to permute
[time cells]=size(df);
if time < cells
    frameShift = randi(time,cells);
else
    frameShift =randperm(time,cells);
end
for cell=1:cells
    shiftDf(:,cell) = circshift(df(:,cell),frameShift(cell),1);
end

