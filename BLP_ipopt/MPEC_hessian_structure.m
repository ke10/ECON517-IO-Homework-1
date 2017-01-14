function hessian_structure = MPEC_hessian_structure(auxdata) 

%%%%%%%%%%%%%%%
% Provides the Sparsity Pattern of the Hessian
%%%%%%%%%%%%%%%

HessianPattern=auxdata{18};
HessianPattern_tril= tril(HessianPattern);
hessian_structure=sparse(HessianPattern_tril);
    
end