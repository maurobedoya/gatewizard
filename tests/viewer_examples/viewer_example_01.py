from gatewizard.core.viewer import MolecularViewer

viewer = MolecularViewer()
print(f"MolecularViewer created: {viewer}")
print(f"Structure loaded: {viewer.structure is not None}")
