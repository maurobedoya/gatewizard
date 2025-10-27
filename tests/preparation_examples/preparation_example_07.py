from gatewizard.core.preparation import PreparationManager

analyzer = PreparationManager()
his_states = analyzer.get_available_states("HIS")
print(his_states)
