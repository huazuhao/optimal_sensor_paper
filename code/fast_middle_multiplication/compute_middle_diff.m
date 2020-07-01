
function [newnew,newold,oldnew] = compute_middle_diff(new_obs,selected_obs,wt,w,total_m,n,rank_cut_off_for_integration)

%%
c_new_m_c_new_w = zeros(n+2,n+2);
c_new_m_c_new_w(new_obs,:) = total_m(new_obs,new_obs)*w(new_obs,:);
first_term = zeros(n+2,n+2);
for row = n+3-rank_cut_off_for_integration:n+2
    for column = n+3-rank_cut_off_for_integration:n+2 
        first_term(row,column) = wt(row,:)*c_new_m_c_new_w(:,column);
    end
end

newnew = first_term;

%%
newold = [];
for column = new_obs-1:new_obs+1
    if ismember(column,selected_obs) == 1
        data_column = [new_obs;column;total_m(new_obs,column)];
        newold = [newold,data_column];
    end
end
c_new_m_c_old_w = zeros(n+2,n+2);
newold_data_size = size(newold);
for data_column = 1:newold_data_size(2)
    data = newold(:,data_column);
    row = data(1);
    column = data(2);
    factor = data(3);
    c_new_m_c_old_w(row,:) = c_new_m_c_old_w(row,:)+factor*w(column,:);
end
second_term = zeros(n+2,n+2);
for row = n+3-rank_cut_off_for_integration:n+2
    for column = n+3-rank_cut_off_for_integration:n+2 
        second_term(row,column) = wt(row,:)*c_new_m_c_old_w(:,column);
    end
end
newold = second_term;

%%
oldnew = [];
for row = new_obs-1:new_obs+1
    if ismember(row,selected_obs)==1
        data_column = [row;new_obs;total_m(row,new_obs)];
        oldnew = [oldnew,data_column];
    end
end

c_old_m_c_new_w = zeros(n+2,n+2);
oldnew_data_size = size(oldnew);
for data_column = 1:oldnew_data_size(2)
    data = oldnew(:,data_column);
    row = data(1);
    column = data(2);
    factor = data(3);
    c_old_m_c_new_w(row,:) = factor*w(column,:); %selecting row of w
end

third_term = zeros(n+2,n+2);
for row = n+3-rank_cut_off_for_integration:n+2
    for column = n+3-rank_cut_off_for_integration:n+2 
        third_term(row,column) = wt(row,:)*c_old_m_c_new_w(:,column);
    end
end
oldnew = third_term;
