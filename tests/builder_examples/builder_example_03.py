from gatewizard.tools.force_fields import ForceFieldManager

ff_manager = ForceFieldManager()

# Get available water models
water_models = ff_manager.get_water_models()

print("\nAvailable Water Models:")
for water in water_models:
    print(f"  - {water}")

print(f"\nTotal: {len(water_models)} water models available")
