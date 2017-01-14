function dc_structure = MPEC_jacobian_structure(auxdata) 

%%%%%%%%%%%%%%%
% Provides the Sparsity Pattern of the Jacobian
%%%%%%%%%%%%%%%

ConsPattern=auxdata{1};

dc_structure=sparse(ConsPattern);
    
end