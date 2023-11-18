% save two variables within parfor loop

function parsave_v2(fname,x,y)
    save(fname,'x','y','-v7.3')
end