function mov_mean = moving_mean(old_data,new_data,cnt)

if cnt == 1
    mov_mean = new_data;
else
    mov_mean = (((cnt-1).*old_data) + new_data)./cnt;
end