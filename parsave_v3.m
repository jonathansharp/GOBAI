% save three variables within parfor loop

function parsave_v3(fname,x,y,z)
    save(fname,'x','y','z','-v7.3')
end