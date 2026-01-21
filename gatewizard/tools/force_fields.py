# gatewizard/tools/force_fields.py
# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Constanza González and Mauricio Bedoya

"""
Force field management and validation for molecular dynamics simulations.

This module provides utilities for managing force field parameters,
validating combinations, and providing recommendations.

Validation References:
[1] Tian, C.; Kasavajhala, K.; Belfon, K. A. A.; Raguette, L.; Huang, H.; Migues, A. N.;
    Bickel, J.; Wang, Y.; Pincay, J.; Wu, Q.; Simmerling, C. ff19SB: Amino-Acid-Specific
    Protein Backbone Parameters Trained against Quantum Mechanics Energy Surfaces in Solution.
    J. Chem. Theory Comput. 2020, 16 (1), 528–552. DOI: 10.1021/acs.jctc.9b00591

[2] Maier, J. A.; Martinez, C.; Kasavajhala, K.; Wickstrom, L.; Hauser, K. E.; Simmerling, C.
    ff14SB: Improving the Accuracy of Protein Side Chain and Backbone Parameters from ff99SB.
    J. Chem. Theory Comput. 2015, 11 (8), 3696–3713. DOI: 10.1021/acs.jctc.5b00255

[3] Debiec, K. T.; Cerutti, D. S.; Baker, L. R.; Gronenborn, A. M.; Case, D. A.; Chong, L. T.
    Further along the Road Less Traveled: AMBER ff15ipq, an Original Protein Force Field Built
    on a Self-Consistent Physical Model. J. Chem. Theory Comput. 2016, 12 (8), 3926–3947.
    DOI: 10.1021/acs.jctc.6b00567

[4] Dickson, C. J.; Walker, R. C.; Gould, I. R. Lipid21: Complex Lipid Membrane Simulations
    with AMBER. J. Chem. Theory Comput. 2022, 18 (3), 1726–1736.
    DOI: 10.1021/acs.jctc.1c01217

[5] Gould, I. R.; Skjevik, A. A.; Dickson, C. J.; Madej, B. D.; Walker, R. C. Lipid17:
    A Comprehensive AMBER Force Field for the Simulation of Zwitterionic and Anionic Lipids.
    Manuscript in preparation, 2018. (This manuscript was never published. Lipid17 parameters are included in Amber since version 17.)

[6] Case, D. A. et al. Amber 2024 Reference Manual; University of California: San Francisco,
    2024. https://ambermd.org/doc12/Amber24.pdf

Validation Status:
- ff14SB + Lipid21 + TIP3P: Fully validated [2,4]
- ff14SB + Lipid17 + TIP3P: Validated (packmol-memgen default) [2,5]
- ff19SB + Lipid21 + OPC: Amber manual recommendation [1,6]
- ff19SB + Lipid21 + TIP3P: Suboptimal for protein [1,4]
- ff15ipq: Only validated for proteins with SPC/Eb [3], no lipid validation
"""

from typing import Dict, List, Set, Tuple, Optional, Any
from dataclasses import dataclass

from gatewizard.utils.logger import get_logger

logger = get_logger(__name__)

@dataclass
class ForceFieldInfo:
    """Information about a force field."""
    name: str
    description: str
    version: str
    compatible_with: List[str]
    recommended_for: List[str]
    year: int
    reference: str = ""
    notes: str = ""

class ForceFieldManager:
    """
    Manager for force field parameters and compatibility.
    
    This class provides information about available force fields,
    their compatibility, and recommendations for different systems.
    """
    
    def __init__(self):
        """Initialize the force field manager."""
        self._water_models = self._init_water_models()
        self._protein_force_fields = self._init_protein_force_fields()
        self._lipid_force_fields = self._init_lipid_force_fields()
        self._available_ions = self._init_available_ions()
        self._available_lipids = self._init_available_lipids()
    
    def _init_water_models(self) -> Dict[str, ForceFieldInfo]:
        """Initialize water model definitions."""
        return {
            "tip3p": ForceFieldInfo(
                name="TIP3P",
                description="Three-point transferable intermolecular potential",
                version="1983",
                # Validated: ff14SB+lipid17 [2,5], ff14SB+lipid21 [2,4]
                # Suboptimal: ff19SB+lipid17/lipid21 (prefer OPC) [1]
                # Not validated: ff15ipq with any lipid FF
                compatible_with=["ff14SB", "lipid17", "lipid21"],
                recommended_for=["general", "membrane"],
                year=1983,
                reference="L. Jorgensen, Jayaraman Chandrasekhar, Jeffry D. Madura, Roger W. Impey, Michael L. Klein; Comparison of simple potential functions for simulating liquid water. J. Chem. Phys. 15 July 1983; 79 (2): 926-935. https://doi.org/10.1063/1.445869",
                notes="Validated for ff14SB with both lipid17 and lipid21. Suboptimal for ff19SB (prefer OPC)."
            ),
            "tip4p": ForceFieldInfo(
                name="TIP4P",
                description="Four-point transferable intermolecular potential",
                version="1983",
                # Not validated with any lipid FF for membrane simulations
                # Protein-only: ff14SB, ff19SB
                compatible_with=["ff14SB", "ff19SB"],
                recommended_for=["general"],
                year=1983,
                reference="L. Jorgensen, Jayaraman Chandrasekhar, Jeffry D. Madura, Roger W. Impey, Michael L. Klein; Comparison of simple potential functions for simulating liquid water. J. Chem. Phys. 15 July 1983; 79 (2): 926-935. https://doi.org/10.1063/1.445869",
                notes="Not validated with lipid force fields. Use for protein-only simulations."
            ),
            "tip4pew": ForceFieldInfo(
                name="TIP4P-Ew",
                description="TIP4P optimized for Ewald summation",
                version="2004",
                # Not validated with any lipid FF for membrane simulations
                # Protein-only: ff14SB, ff19SB
                compatible_with=["ff14SB", "ff19SB"],
                recommended_for=["PME"],
                year=2004,
                reference="Horn, H. W., Swope, W. C., Pitera, J. W., Madura, J. D., Dick, T. J., Hura, G. L., & Head-Gordon, T. (2004). Development of an improved four-site water model for biomolecular simulations: TIP4P-Ew. The Journal of chemical physics, 120(20), 9665-9678. https://doi.org/10.1063/1.1683075",
                notes="Not validated with lipid force fields. Use for protein-only simulations."
            ),
            "spce": ForceFieldInfo(
                name="SPC/E",
                description="Extended simple point charge model",
                version="1987",
                # Not validated with any lipid FF for membrane simulations
                # Protein-only: ff14SB, ff19SB
                compatible_with=["ff14SB", "ff19SB"],
                recommended_for=["general"],
                year=1987,
                reference="Berendsen, H. J. C., Grigera, J. R., & Straatsma, T. P. (1987). The missing term in effective pair potentials. Journal of Physical Chemistry, 91(24), 6269-6271. https://doi.org/10.1021/j100308a038",
                notes="Not validated with lipid force fields. Use for protein-only simulations."
            ),
            "opc": ForceFieldInfo(
                name="OPC",
                description="Optimal point charge water model",
                version="2014",
                # Amber manual recommendation: ff19SB+lipid21 [1,6]
                # Not validated: ff14SB+lipid17/lipid21, ff19SB+lipid17, ff15ipq+any lipid
                compatible_with=["ff19SB", "lipid21"],
                recommended_for=["protein", "latest", "ff19SB"],
                year=2014,
                reference="Izadi, S., Anandakrishnan, R., & Onufriev, A.V. (2014). Building Water Models: A Different Approach. The Journal of Physical Chemistry Letters, 5, 3863 - 3871. https://pubs.acs.org/doi/10.1021/jz501780a",
                notes="Amber manual recommends OPC with ff19SB+lipid21 [6]. Not validated with ff14SB or lipid17."
            ),
            "opc3": ForceFieldInfo(
                name="OPC3",
                description="Three-point optimal point charge model",
                version="2016",
                # Not validated with any lipid FF for membrane simulations
                # Protein-only: ff14SB, ff19SB
                compatible_with=["ff14SB", "ff19SB"],
                recommended_for=["protein"],
                year=2016,
                reference="Saeed Izadi, Alexey V. Onufriev; Accuracy limit of rigid 3-point water models. J. Chem. Phys. 21 August 2016; 145 (7): 074501. https://doi.org/10.1063/1.4960175",
                notes="Not validated with lipid force fields. Use for protein-only simulations."
            ),
            "spceb": ForceFieldInfo(
                name="SPC/Eb",
                description="SPC/E optimized for biomolecules",
                version="2010",
                # Only validated with ff15ipq for protein simulations [3]
                # Not validated with lipid FFs
                compatible_with=["ff15ipq"],
                recommended_for=["protein", "ff15ipq"],
                year=2010,
                reference="Takemura, K., & Kitao, A. (2012). Water model tuning for improved reproduction of rotational diffusion and NMR spectral density. The journal of physical chemistry. B, 116(22), 6279-6287. https://doi.org/10.1021/jp301100g",
                notes="Specifically validated with ff15ipq for protein simulations [3]. Not validated with lipid force fields."
            ),
            "fb3": ForceFieldInfo(
                name="FB3",
                description="Force-balanced three-point model",
                version="2019",
                compatible_with=["ff19SB"],
                recommended_for=["protein"],
                year=2019,
                reference="Wang, L. P., Martinez, T. J., & Pande, V. S. (2014). Building Force Fields: An Automatic, Systematic, and Reproducible Approach. The journal of physical chemistry letters, 5(11), 1885-1891. https://doi.org/10.1021/jz500737m"
            )
        }
    
    def _init_protein_force_fields(self) -> Dict[str, ForceFieldInfo]:
        """Initialize protein force field definitions."""
        return {
            "ff14SB": ForceFieldInfo(
                name="ff14SB",
                description="Amber force field with improved side-chain torsions",
                version="2014",
                # Validated with lipids: tip3p+lipid17 [2,5], tip3p+lipid21 [2,4]
                # Protein-only: tip4p, tip4pew, spce, opc, opc3
                compatible_with=["tip3p", "tip4p", "tip4pew", "spce", "opc", "opc3", "lipid17", "lipid21"],
                recommended_for=["protein", "general", "membrane"],
                year=2014,
                reference="Maier, J. A.; Martinez, C.; Kasavajhala, K.; Wickstrom, L.; Hauser, K. E.; Simmerling, C. ff14SB: Improving the Accuracy of Protein Side Chain and Backbone Parameters from ff99SB. J. Chem. Theory Comput. 2015, 11 (8), 3696–3713. DOI: 10.1021/acs.jctc.5b00255",
                notes="Fully validated for membrane simulations with TIP3P water and lipid17/lipid21 [2,4,5]. Use TIP3P for membrane systems."
            ),
            "ff15ipq": ForceFieldInfo(
                name="ff15ipq",
                description="Improved protein backbone parameters",
                version="2015",
                # Only validated with SPC/Eb for protein simulations [3]
                # NOT validated with any lipid force field
                compatible_with=["spceb"],
                recommended_for=["protein", "folding"],
                year=2015,
                reference="Debiec, K. T.; Cerutti, D. S.; Baker, L. R.; Gronenborn, A. M.; Case, D. A.; Chong, L. T. Further along the Road Less Traveled: AMBER ff15ipq, an Original Protein Force Field Built on a Self-Consistent Physical Model. J. Chem. Theory Comput. 2016, 12 (8), 3926–3947. DOI: 10.1021/acs.jctc.6b00567",
                notes="Only validated with SPC/Eb water for protein simulations [3]. NOT validated with any lipid force field - do not use for membrane simulations."
            ),
            "ff19SB": ForceFieldInfo(
                name="ff19SB",
                description="Latest Amber protein force field with optimized backbone",
                version="2019",
                # Amber manual recommends: opc+lipid21 [1,6]
                # Suboptimal: tip3p (prefer OPC) [1]
                # Protein-only: tip4p, tip4pew, spce, opc3, fb3
                compatible_with=["tip3p", "tip4p", "tip4pew", "spce", "opc", "opc3", "fb3", "lipid17", "lipid21"],
                recommended_for=["protein", "latest"],
                year=2019,
                reference="Tian, C.; Kasavajhala, K.; Belfon, K. A. A.; Raguette, L.; Huang, H.; Migues, A. N.; Bickel, J.; Wang, Y.; Pincay, J.; Wu, Q.; Simmerling, C. ff19SB: Amino-Acid-Specific Protein Backbone Parameters Trained against Quantum Mechanics Energy Surfaces in Solution. J. Chem. Theory Comput. 2020, 16 (1), 528–552. DOI: 10.1021/acs.jctc.9b00591",
                notes="Amber manual recommends OPC water with lipid21 for membrane simulations [1,6]. TIP3P is suboptimal for ff19SB."
            )
        }
    
    def _init_lipid_force_fields(self) -> Dict[str, ForceFieldInfo]:
        """Initialize lipid force field definitions."""
        return {
            "lipid17": ForceFieldInfo(
                name="Lipid17",
                description="Amber lipid force field compatible with protein force fields",
                version="2017",
                # Validated: ff14SB+tip3p [2,5] (packmol-memgen default)
                # Suboptimal: ff19SB+tip3p (prefer lipid21+OPC) [1]
                # Not validated: ff15ipq with any water, ff19SB+OPC
                compatible_with=["ff14SB", "ff19SB", "tip3p"],
                recommended_for=["membrane", "lipid"],
                year=2017,
                reference="There was a manuscript in preparation by Gould et al. in 2018, but it was never published. Lipid17 parameters are included in Amber since version 17.",
                notes="Validated with ff14SB+TIP3P (packmol-memgen default) [2,5]. For ff19SB, prefer lipid21."
            ),
            "lipid21": ForceFieldInfo(
                name="Lipid21",
                description="Updated Amber lipid force field with improved parameters",
                version="2021",
                # Validated: ff14SB+tip3p [2,4]
                # Recommended: ff19SB+opc [1,6] (Amber manual)
                # Suboptimal: ff19SB+tip3p [1]
                # Not validated: ff15ipq with any water
                compatible_with=["ff14SB", "ff19SB", "tip3p", "opc"],
                recommended_for=["membrane", "lipid", "latest"],
                year=2021,
                reference="Dickson, C. J.; Walker, R. C.; Gould, I. R. Lipid21: Complex Lipid Membrane Simulations with AMBER. J. Chem. Theory Comput. 2022, 18 (3), 1726–1736. DOI: 10.1021/acs.jctc.1c01217",
                notes="Validated with ff14SB+TIP3P [2,4]. Amber manual recommends ff19SB+OPC for latest combination [1,6]."
            )
        }
    
    def _init_available_ions(self) -> Dict[str, Dict[str, List[Tuple[str, int]]]]:
        """Initialize available ions with their charges."""
        return {
            'cations': {
                'monovalent': [
                    ('Li+', 1), ('Na+', 1), ('K+', 1), ('Rb+', 1), ('Cs+', 1),
                    ('Ag+', 1), ('Cu+', 1), ('Tl+', 1), ('H3O+', 1), ('NH4+', 1)
                ],
                'divalent': [
                    ('Be2+', 2), ('Mg2+', 2), ('Ca2+', 2), ('Sr2+', 2), ('Ba2+', 2),
                    ('Ra2+', 2), ('Mn2+', 2), ('Fe2+', 2), ('Co2+', 2), ('Ni2+', 2),
                    ('Cu2+', 2), ('Zn2+', 2), ('Cd2+', 2), ('Hg2+', 2), ('Pb2+', 2),
                    ('Sn2+', 2), ('Pd2+', 2), ('Pt2+', 2), ('V2+', 2), ('Cr2+', 2),
                    ('Eu2+', 2), ('Sm2+', 2), ('Yb2+', 2)
                ],
                'trivalent': [
                    ('Al3+', 3), ('Cr3+', 3), ('Fe3+', 3), ('In3+', 3), ('Tl3+', 3),
                    ('Y3+', 3), ('La3+', 3), ('Ce3+', 3), ('Pr3+', 3), ('Nd3+', 3),
                    ('Sm3+', 3), ('Eu3+', 3), ('Gd3+', 3), ('Tb3+', 3), ('Dy3+', 3),
                    ('Er3+', 3), ('Tm3+', 3), ('Lu3+', 3)
                ],
                'tetravalent': [
                    ('Zr4+', 4), ('Hf4+', 4), ('Ce4+', 4), ('Th4+', 4), ('U4+', 4), ('Pu4+', 4)
                ]
            },
            'anions': {
                'monovalent': [
                    ('F-', -1), ('Cl-', -1), ('Br-', -1), ('I-', -1)
                ]
            }
        }
    
    def _init_available_lipids(self) -> List[str]:
        """Initialize list of available lipids."""
        return [
            # Phosphatidylcholine (PC) lipids
            "DHPC", "DLPC", "DMPC", "DPPC", "DSPC", "DOPC", "POPC", "PLPC", "SOPC",
            
            # Phosphatidylethanolamine (PE) lipids
            "DHPE", "DLPE", "DMPE", "DPPE", "DSPE", "DOPE", "POPE", "PLPE", "SOPE",
            
            # Phosphatidylserine (PS) lipids
            "DHPS", "DLPS", "DMPS", "DPPS", "DSPS", "DOPS", "POPS", "PLPS", "SOPS",
            
            # Phosphatidylglycerol (PG) lipids
            "DHPG", "DLPG", "DMPG", "DPPG", "DSPG", "DOPG", "POPG", "PLPG", "SOPG",
            
            # Phosphatidic acid (PA) lipids
            "DHPA", "DLPA", "DMPA", "DPPA", "DSPA", "DOPA", "POPA", "PLPA", "SOPA",
            
            # Sphingomyelins and other lipids
            "PSM", "ASM", "LSM", "MSM", "HSM", "OSM",
            
            # Cholesterol
            "CHL1",
            
            # Other specialized lipids
            "CARDIOLIPIN", "PIP", "PIP2"
        ]
    
    def get_water_models(self) -> List[str]:
        """Get list of available water models."""
        return list(self._water_models.keys())
    
    def get_protein_force_fields(self) -> List[str]:
        """Get list of available protein force fields."""
        return list(self._protein_force_fields.keys())
    
    def get_lipid_force_fields(self) -> List[str]:
        """Get list of available lipid force fields."""
        return list(self._lipid_force_fields.keys())
    
    def get_available_lipids(self) -> List[str]:
        """Get list of available lipids."""
        return sorted(self._available_lipids)
    
    def get_available_cations(self) -> List[str]:
        """Get list of available cations."""
        cations = []
        for category in self._available_ions['cations'].values():
            cations.extend([ion[0] for ion in category])
        return sorted(cations)
    
    def get_available_anions(self) -> List[str]:
        """Get list of available anions."""
        anions = []
        for category in self._available_ions['anions'].values():
            anions.extend([ion[0] for ion in category])
        return sorted(anions)
    
    def validate_combination(
        self, 
        water_model: str, 
        protein_ff: str, 
        lipid_ff: str
    ) -> Tuple[bool, str, bool]:
        """
        Validate a force field combination.
        
        Args:
            water_model: Water model name
            protein_ff: Protein force field name
            lipid_ff: Lipid force field name
            
        Returns:
            Tuple of (is_valid, message, is_warning)
            - is_valid: False only for unknown components, True otherwise
            - message: Validation message
            - is_warning: True if combination is unvalidated but allowed
        """
        # Check if components exist (these are actual errors)
        if water_model not in self._water_models:
            return False, f"Unknown water model: {water_model}", False
        
        if protein_ff not in self._protein_force_fields:
            return False, f"Unknown protein force field: {protein_ff}", False
        
        if lipid_ff not in self._lipid_force_fields:
            return False, f"Unknown lipid force field: {lipid_ff}", False
        
        # Check compatibility (these are warnings, not errors)
        water_info = self._water_models[water_model]
        protein_info = self._protein_force_fields[protein_ff]
        lipid_info = self._lipid_force_fields[lipid_ff]
        
        warnings = []
        
        # Check water-protein compatibility
        if protein_ff not in water_info.compatible_with:
            warnings.append(f"Water model {water_model} and protein force field {protein_ff}")
        
        # Check water-lipid compatibility
        if lipid_ff not in water_info.compatible_with:
            warnings.append(f"Water model {water_model} and lipid force field {lipid_ff}")
        
        # Check protein-lipid compatibility
        if lipid_ff not in protein_info.compatible_with:
            warnings.append(f"Protein force field {protein_ff} and lipid force field {lipid_ff}")
        
        if warnings:
            warning_msg = (
                "WARNING: The following combination(s) have not been validated in the literature:\n\n"
                + "\n".join(f"  - {w}" for w in warnings)
                + "\n\nThis combination may produce unreliable results. "
                + "You may proceed at your own risk for testing/development purposes."
            )
            return True, warning_msg, True
        
        return True, "Force field combination is valid and well-tested", False
    
    def get_recommendations(self, system_type: str = "membrane") -> Dict[str, str]:
        """
        Get recommended force field combination for a system type.
        
        Args:
            system_type: Type of system (membrane, protein, general, latest)
            
        Returns:
            Dictionary with recommended force fields
        """
        if system_type == "membrane":
            # Fully validated combination [2,4]
            return {
                "water_model": "tip3p",
                "protein_ff": "ff14SB",
                "lipid_ff": "lipid21",
                "reason": "Fully validated combination for membrane protein simulations [Refs 2,4]"
            }
        elif system_type == "protein":
            return {
                "water_model": "opc",
                "protein_ff": "ff19SB",
                "lipid_ff": "lipid21",
                "reason": "Latest protein force field with recommended water model [Refs 1,6]"
            }
        elif system_type == "latest":
            # Amber manual recommendation [1,6]
            return {
                "water_model": "opc",
                "protein_ff": "ff19SB",
                "lipid_ff": "lipid21",
                "reason": "Amber manual recommended combination for modern simulations [Refs 1,6]"
            }
        else:  # general
            # packmol-memgen default [2,5]
            return {
                "water_model": "tip3p",
                "protein_ff": "ff14SB", 
                "lipid_ff": "lipid17",
                "reason": "Validated general-purpose combination (packmol-memgen default) [Refs 2,5]"
            }
    
    def get_force_field_info(self, ff_name: str, ff_type: str) -> Optional[ForceFieldInfo]:
        """
        Get detailed information about a force field.
        
        Args:
            ff_name: Force field name
            ff_type: Type (water, protein, lipid)
            
        Returns:
            ForceFieldInfo object or None if not found
        """
        if ff_type == "water":
            return self._water_models.get(ff_name)
        elif ff_type == "protein":
            return self._protein_force_fields.get(ff_name)
        elif ff_type == "lipid":
            return self._lipid_force_fields.get(ff_name)
        else:
            return None
    
    def validate_lipid(self, lipid_name: str) -> bool:
        """
        Validate if a lipid is available.
        
        Args:
            lipid_name: Lipid name to validate
            
        Returns:
            True if lipid is available, False otherwise
        """
        return lipid_name in self._available_lipids
    
    def validate_ion(self, ion_name: str) -> Tuple[bool, int]:
        """
        Validate an ion and return its charge.
        
        Args:
            ion_name: Ion name to validate
            
        Returns:
            Tuple of (is_valid, charge)
        """
        # Check cations
        for category in self._available_ions['cations'].values():
            for ion, charge in category:
                if ion == ion_name:
                    return True, charge
        
        # Check anions
        for category in self._available_ions['anions'].values():
            for ion, charge in category:
                if ion == ion_name:
                    return True, charge
        
        return False, 0