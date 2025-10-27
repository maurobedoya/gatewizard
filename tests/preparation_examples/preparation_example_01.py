from gatewizard.core.preparation import PreparationManager

analyzer = PreparationManager(propka_version="3")
print(f"Using PROPKA version: {analyzer.propka_version}")
