function struct_out = extractDay1(struct_in)

temp = fields(struct_in);
for i=1:length(temp)
    struct_out.(temp{i}) = struct_in.(temp{i})(:,1);
end

end