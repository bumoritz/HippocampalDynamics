function [struct_out] = createStimSquare(data_in,d_info,shape)

% initialise struct_out
if strcmp(shape,'seq-ctrl_switch-postswitch')
    struct_out.seq_switch = [];
    struct_out.seq_postswitch = [];
    struct_out.ctrl_switch = [];
    struct_out.ctrl_postswitch = [];
end

if strcmp(shape,'seq-ctrl_switch-postswitch')
    for i=1:d_info.numAnimals
        if d_info.group(i)==7
            struct_out.seq_switch = [struct_out.seq_switch; data_in(i,2)];
            struct_out.seq_postswitch = [struct_out.seq_postswitch; data_in(i,3)];
            struct_out.ctrl_switch = [struct_out.ctrl_switch; data_in(i,4)];
            struct_out.ctrl_postswitch = [struct_out.ctrl_postswitch; data_in(i,5)];
        elseif d_info.group(i)==8
            struct_out.seq_switch = [struct_out.seq_switch; data_in(i,4)];
            struct_out.seq_postswitch = [struct_out.seq_postswitch; data_in(i,5)];
            struct_out.ctrl_switch = [struct_out.ctrl_switch; data_in(i,2)];
            struct_out.ctrl_postswitch = [struct_out.ctrl_postswitch; data_in(i,3)];
        end
    end
end

end

