function []= testingGPU ()
tTotal = tic;
    for  i=1:1000000
        gpuArray.rand(1000,1,'single');
    end
    toc(tTotal);
    
    
    %tTotal = tic;    
    %for  i=1:10000
    %    rand(1000000,1,'single');
    %end
        
    
%toc(tTotal);
end
