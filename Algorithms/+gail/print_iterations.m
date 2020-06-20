function print_iterations(irep, desc, oneline)
howoftenrep = 10^floor(log10(irep));
if mod(irep, howoftenrep)==0
  if (oneline==true)
    if irep == 1
      fprintf('\n%s = %d', desc, irep);
    else
      fprintf(' ... %d', irep);
    end
  else
    fprintf('irep = %d\n', irep);
  end
end
