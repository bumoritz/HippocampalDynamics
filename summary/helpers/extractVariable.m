function extracted = extractVariable(struct_in,variable,output_type,entries)
% struct_in = d; variable = 'paq_beh.numFrameTriggers_original'; output_type = 'array'; entries = 'first';
% struct_in = d; variable = 'paq_beh.licks_od2'; output_type = 'cell'; entries = 'firstOfEach';
% struct_in = paq_beh; variable = 'licks_od2'; output_type = 'array'; entries = 'number';
% output_types: array, cell; entries: first, all, firstOfEach, number (implement mean/min/max) 

variable_parts = split(variable,'.');
variable_numParts = length(variable_parts);

if strcmp(output_type,'array')
    extracted = nan(size(struct_in));
elseif strcmp(output_type,'cell')
	extracted = {};
end
for i=1:size(struct_in,1)
    for j=1:size(struct_in,2)
        
        move_on = false;
        if ~(size(struct_in,1)==1 && size(struct_in,2)==1)
            if variable_numParts==1
                move_on = isfield(struct_in{i,j},variable_parts{1});
            elseif variable_numParts==2
                move_on = isfield(struct_in{i,j},variable_parts{1}) && isfield(struct_in{i,j}.(variable_parts{1}),variable_parts{2});
            elseif variable_numParts==3
                move_on = (isfield(struct_in{i,j},variable_parts{1}) && isfield(struct_in{i,j}.(variable_parts{1}),variable_parts{2})) && isfield(struct_in{i,j}.(variable_parts{1}).(variable_parts{2}),variable_parts{3});
            elseif variable_numParts==4
                move_on = ((isfield(struct_in{i,j},variable_parts{1}) && isfield(struct_in{i,j}.(variable_parts{1}),variable_parts{2})) && isfield(struct_in{i,j}.(variable_parts{1}).(variable_parts{2}),variable_parts{3})) && isfield(struct_in{i,j}.(variable_parts{1}).(variable_parts{2}).(variable_parts{3}),variable_parts{4});
            elseif variable_numParts==5
                move_on = (((isfield(struct_in{i,j},variable_parts{1}) && isfield(struct_in{i,j}.(variable_parts{1}),variable_parts{2})) && isfield(struct_in{i,j}.(variable_parts{1}).(variable_parts{2}),variable_parts{3})) && isfield(struct_in{i,j}.(variable_parts{1}).(variable_parts{2}).(variable_parts{3}),variable_parts{4})) && isfield(struct_in{i,j}.(variable_parts{1}).(variable_parts{2}).(variable_parts{3}).(variable_parts{4}),variable_parts{5});
            end
            if variable_numParts==1 && move_on
                content = struct_in{i,j}.(variable_parts{1});
            elseif variable_numParts==2 && move_on
                content = struct_in{i,j}.(variable_parts{1}).(variable_parts{2});
            elseif variable_numParts==3 && move_on
                content = struct_in{i,j}.(variable_parts{1}).(variable_parts{2}).(variable_parts{3});
            elseif variable_numParts==4 && move_on
                content = struct_in{i,j}.(variable_parts{1}).(variable_parts{2}).(variable_parts{3}).(variable_parts{4});
            elseif variable_numParts==5 && move_on
                content = struct_in{i,j}.(variable_parts{1}).(variable_parts{2}).(variable_parts{3}).(variable_parts{4}).(variable_parts{5});
            end
        else
            if variable_numParts==1
                move_on = isfield(struct_in,variable_parts{1});
            elseif variable_numParts==2
                move_on = isfield(struct_in,variable_parts{1}) && isfield(struct_in.(variable_parts{1}),variable_parts{2});
            elseif variable_numParts==3
                move_on = (isfield(struct_in,variable_parts{1}) && isfield(struct_in.(variable_parts{1}),variable_parts{2})) && isfield(struct_in.(variable_parts{1}).(variable_parts{2}),variable_parts{3});
            elseif variable_numParts==4
                move_on = ((isfield(struct_in,variable_parts{1}) && isfield(struct_in.(variable_parts{1}),variable_parts{2})) && isfield(struct_in.(variable_parts{1}).(variable_parts{2}),variable_parts{3})) && isfield(struct_in.(variable_parts{1}).(variable_parts{2}).(variable_parts{3}),variable_parts{4});
            elseif variable_numParts==5
                move_on = (((isfield(struct_in,variable_parts{1}) && isfield(struct_in.(variable_parts{1}),variable_parts{2})) && isfield(struct_in.(variable_parts{1}).(variable_parts{2}),variable_parts{3})) && isfield(struct_in.(variable_parts{1}).(variable_parts{2}).(variable_parts{3}),variable_parts{4})) && isfield(struct_in.(variable_parts{1}).(variable_parts{2}).(variable_parts{3}).(variable_parts{4}),variable_parts{5});
            end
            if variable_numParts==1 && move_on
                content = struct_in.(variable_parts{1});
            elseif variable_numParts==2 && move_on
                content = struct_in.(variable_parts{1}).(variable_parts{2});
            elseif variable_numParts==3 && move_on
                content = struct_in.(variable_parts{1}).(variable_parts{2}).(variable_parts{3});
            elseif variable_numParts==4 && move_on
                content = struct_in.(variable_parts{1}).(variable_parts{2}).(variable_parts{3}).(variable_parts{4});
            elseif variable_numParts==5 && move_on
                content = struct_in.(variable_parts{1}).(variable_parts{2}).(variable_parts{3}).(variable_parts{4}).(variable_parts{5});
            end
        end
        
        
        if move_on

            if strcmp(output_type,'array')
                if strcmp(entries,'first')
                    extracted(i,j) = content(1);
                elseif strcmp(entries,'single') % used in respXperfSummary
                    extracted(i,j) = content;
                elseif strcmp(entries,'firstOfEach') % used in lickingAnalysis
                    extracted = nan(length(content),1);
                    for t=1:length(content)
                        if ~isempty(content{t})
                            extracted(t) = content{t}(1);
                        end
                    end
                elseif strcmp(entries,'number') % used in lickingAnalysis
                    extracted = zeros(length(content),1);
                    for t=1:length(content)
                        if ~isempty(content{t})
                            extracted(t) = length(content{t});
                        end
                    end
                elseif strcmp(entries,'all')
                    extracted = content;
                end

            elseif strcmp(output_type,'cell')
                if strcmp(entries,'first')
                    extracted{i,j} = content(1);
                elseif strcmp(entries,'firstOfEach')
                    extracted{i,j} = nan(length(content),1);
                    for t=1:length(content)
                        if ~isempty(content{t})
                            extracted{i,j}(t) = content{t}(1);
                        end
                    end
                elseif strcmp(entries,'number')
                    extracted{i,j} = zeros(length(content),1);
                    for t=1:length(content)
                        if ~isempty(content{t})
                            extracted{i,j}(t) = length(content{t});
                        end
                    end
                elseif strcmp(entries,'all')
                    extracted{i,j} = content;
                end
            end

        end
    end
end