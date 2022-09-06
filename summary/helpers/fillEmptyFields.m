function struct_out = fillEmptyFields(struct_in,fieldname,templateSession_noswitch,templateSession_switch)
% struct_in = d; fieldname = 'perf';

% define source sessions as templates for creating skaffolds
if nargin < 3
    templateSession_noswitch = [23,5];  % default: Maguire (2nd switch)
    templateSession_switch = [23,4]; % default: Maguire (2nd switch)
end

% create skaffold for no-switch sessions
template_noswitch = struct_in{templateSession_noswitch(1),templateSession_noswitch(2)}.(fieldname);
if strcmp(fieldname,'perf')
    temp = fields(template_noswitch);
    skaffold_noswitch = {};
    for i=1:length(temp)
        temp1 = fields(template_noswitch.(char(temp(i))));
        for j=1:length(temp1)
            temp2 = template_noswitch.(char(temp(i))).(char(temp1(j)));
            temp2(:) = NaN;
            skaffold_noswitch.(fieldname).(char(temp(i))).(char(temp1(j))) = temp2;
        end
    end
end


% create skaffold for switch sessions
template_switch = struct_in{templateSession_switch(1),templateSession_switch(2)}.(fieldname);
if strcmp(fieldname,'perf')
    temp = fields(template_switch);
    skaffold_switch = {};
    for i=1:length(temp)
        temp1 = fields(template_switch.(char(temp(i))));
        for j=1:length(temp1)
            temp2 = template_switch.(char(temp(i))).(char(temp1(j)));
            temp2(:) = NaN;
            skaffold_switch.(fieldname).(char(temp(i))).(char(temp1(j))) = temp2;
        end
    end
end

% paste skaffold into structure
struct_out = struct_in;
for i=1:size(struct_in,1)
    for j=1:size(struct_in,2)
        if isempty(struct_in{i,j})
            if (j==2 | j==4)
                struct_out{i,j} = skaffold_switch;
            else
                struct_out{i,j} = skaffold_noswitch;
            end
        end
    end
end


%% old version

% % create skaffold for no-switch sessions
% template_noswitch = struct_in{templateSession_noswitch(1),templateSession_noswitch(2)};
% skaffold_noswitch = nan(size(template_noswitch));
% 
% % create skaffold for switch sessions
% template_switch = struct_in{templateSession_switch(1),templateSession_switch(2)};
% skaffold_switch = nan(size(template_switch));
% 
% % paste skaffold into structure
% struct_out = struct_in;
% for i=1:size(struct_in,1)
%     for j=1:size(struct_in,2)
%         if isempty(struct_in{i,j})
%             if (j==2 | j==4)
%                 struct_out{i,j} = skaffold_switch;
%             else
%                 struct_out{i,j} = skaffold_noswitch;
%             end
%         end
%     end
% end

end