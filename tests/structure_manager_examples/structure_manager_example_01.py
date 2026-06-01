from gatewizard.core.structure_manager import StructureManager

viewer = StructureManager()
print(f"StructureManager created: {viewer}")
print(f"Structure loaded: {viewer.structure is not None}")
