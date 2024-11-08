clear;clc;close all

obj0 = importdata('his_obj_dps0.csv');
obj1 = importdata('his_obj_dps1.csv');
hobj = importdata('his_obj_actual.csv');

cnt0 = 0;
for i = 1:size(obj0, 1)
    cnt = sum(obj0(i,:) < hobj);
    if cnt == 4
        cnt0 = cnt0 + 1;
    end
end

cnt1 = 0;
for i = 1:size(obj1, 1)
    cnt = sum(obj1(i,:) < hobj);
    if cnt == 4
        cnt1 = cnt1 + 1;
    end
end

% tdata0 = importdata('re0.reference');
% obj0 = tdata0(:,end-3:end);
% tdata1 = importdata('re1.reference');
% obj1 = tdata1(:,end-3:end);

obj0(:, 5) = 0;
obj1(:, 5) = 1;

hobj(:, 5) = 2;
total = [obj0; obj1; hobj];

% total = [obj0; obj1];

l0 = -1;
l1 = size(total, 1);

mark = 0;
while l0 ~= l1
    l0 = size(total, 1);
    mask = zeros(l0, 1);
    for i = 1:l0
        for j = 1:l0
            cnt = sum(total(i, 1:4) < total(j, 1:4));
            if cnt == 4
                mask(j) = mask(j) + 1;
            end
        end
    end
    total(mask~=0, :) = [];
    l1 = size(total, 1);
    mark = mark + 1
end

% ---------------------------------------------------------------
% Case                         | obj0 [%]       | obj1 [%]      |
% ---------------------------------------------------------------
% Synthetic Inflow             | 55.97          | 44.03         |
% ---------------------------------------------------------------
% Historical Inflow (w/ hobj)  | 42.75          | 56.49         |
% ---------------------------------------------------------------
% Historical Inflow (w/ hobj)  | 43.08          | 56.92         |
% ---------------------------------------------------------------
% H0: perfect and binary policies are equally good
% p = 0.02 for merged synthetic results to not reject H0
% p = 0.07 for historically validated results to not reject H0
% based on 1e6 MC experiments
% Only 1 and 2 policies are dominated by historical operations when
% validated using historical inflow from obj0 and obj1, respectively.
% None of both policies dominates historical operations. 