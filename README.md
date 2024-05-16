Analysis Code:
1) ROC analysis
   - simulateROC_Wrapper: Wrapper to run simulateROC with customizable params
   - simulateROC : function to run ROC analysis
2) Shared neural space via PLSC
   - plsc: function to perform PLSC transformation
   - computeSharedNullDistribution: compute null distribution for PLSC analysis and identify number of significant dimensions
   - getNeuralSpace: retrieve shared and unique neural space
3) Coordinated behavior space via CCA
   - computeCoordNullDistribution: compute null distribution for CCA analysis and identify number of significant dimensions
   - getBehaviorSpace: retrieve coordinated and uncoordinated behavior space
4) Compute non-redundant variance explained using Partial Least Square Regression
   - computeNonRedundantVar: compute the non-redundant contribution of specific groups of variables of interest in a model using PLSR and permutation approach
