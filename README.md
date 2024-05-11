Analysis Code for Zhang-Phi 2024
1) ROC analysis
   - simulateROC_Wrapper: Wrapper to run simulateROC with customizable params
   - simulateROC : function to run ROC analysis
2) Decoder analysis
   - to be updated
3) Shared neural space via PLSC
   - plsc: function to perform PLSC transformation
   - computeSharedNullDistribution: compute null distribution for PLSC analysis and identify number of significant dimensions
   - getNeuralSpace: retrieve shared and unique neural space
4) Coordinated behavior space via PLSC
   - computeCoordNullDistribution: compute null distribution for CCA analysis and identify number of significant dimensions
   - getBehaviorSpace: retrieve coordinated and uncoordinated behavior space
5) Non-redundant contribution via PLSR
   - to be updated
