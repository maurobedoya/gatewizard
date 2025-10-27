#!/bin/bash
    # Run preparation script for Gatewizard
    # Generated: 2025-10-19T18:40:16.547815
    
    # New Gatewizard workflow:
    # Stage 1: Pack protein into membrane WITH --parametrize (packmol-memgen) for CRYST1 info
    # Stage 2: Run pdb4amber to prepare PDB file for Amber compatibility  
    # Stage 3: Run tleap to parametrize the system using prepared PDB
    #
    # The bilayer*_lipid.pdb file from Stage 1 provides essential CRYST1 box information
    # This workflow produces more reliable .prmtop/.inpcrd files for MD simulations

    # Save PID for tracking
    echo $$ > "systems/popc_membrane/process.pid"

    # Change to job directory
    cd "systems/popc_membrane"

    # Create logs directory
    mkdir -p logs

    # Define steps for monitoring - dynamic based on configuration
    preoriented=true
    parametrize_requested=true
    
    if [ "$preoriented" = "true" ]; then
        if [ "$parametrize_requested" = "true" ]; then
            declare -a STEPS=("Packmol" "pdb4amber" "tleap")
        else
            declare -a STEPS=("Packmol")
        fi
    else
        if [ "$parametrize_requested" = "true" ]; then
            declare -a STEPS=("MEMEMBED" "Packmol" "pdb4amber" "tleap")
        else
            declare -a STEPS=("MEMEMBED" "Packmol")
        fi
    fi
    
    echo "INFO: Workflow steps: ${STEPS[@]}" | tee -a logs/preparation.log

    # Create Python status utility script
    cat > status_utils.py << 'EOF'
import json
import sys
import os
from datetime import datetime

def update_status(step, msg):
    try:
        # Debug: Print current working directory
        print('DEBUG: Working in directory:', os.getcwd())
        
        # Read current status
        try:
            with open('status.json', 'r') as f:
                status = json.load(f)
            print('DEBUG: Read existing status file')
        except Exception as e:
            print('DEBUG: Creating new status file, error was:', str(e))
            status = {
                'status': 'running',
                'start_time': datetime.now().isoformat(),
                'current_step': 0,
                'total_steps': 5,
                'steps_completed': [],
                'last_update': None
            }

        # Update status
        status['current_step'] = step
        status['last_update'] = datetime.now().isoformat()
        
        # Update message for current step
        if step <= len(status.get('step_messages', [])):
            if 'step_messages' not in status:
                status['step_messages'] = []
            while len(status['step_messages']) < step:
                status['step_messages'].append("")
            status['step_messages'][step-1] = msg
        
        # Write updated status
        with open('status.json', 'w') as f:
            json.dump(status, f, indent=2)
        
        print(f'DEBUG: Updated status - Step {step}: {msg}')
        
    except Exception as e:
        print('ERROR updating status:', str(e))
        import traceback
        traceback.print_exc()

def mark_complete():
    try:
        print('DEBUG: Marking job complete in', os.getcwd())
        
        with open('status.json', 'r') as f:
            status = json.load(f)

        status['status'] = 'completed'
        status['end_time'] = datetime.now().isoformat()
        status['current_step'] = max(status.get('current_step', 0), 5)  # Ensure final step
        
        # Add completion to steps if not already there
        if 'Completed' not in status.get('steps_completed', []):
            status['steps_completed'].append('Completed')

        with open('status.json', 'w') as f:
            json.dump(status, f, indent=2)
            
        print('DEBUG: Job marked as complete')
    except Exception as e:
        print('ERROR marking complete:', str(e))

def handle_error(error_msg):
    try:
        with open('status.json', 'r') as f:
            status = json.load(f)

        status['status'] = 'error'
        status['error'] = error_msg
        status['end_time'] = datetime.now().isoformat()

        with open('status.json', 'w') as f:
            json.dump(status, f, indent=2)
            
        print('ERROR: Job marked as failed -', error_msg)
    except Exception as e:
        print('ERROR handling error:', str(e))

if __name__ == "__main__":
    action = sys.argv[1]
    if action == "update":
        update_status(int(sys.argv[2]), sys.argv[3])
    elif action == "complete":
        mark_complete()
    elif action == "error":
        handle_error(sys.argv[2])
EOF

    # Function to update status
    update_status() {
        step=$1
        msg=$2
        python3 status_utils.py update "$step" "$msg" 2>&1 | tee -a logs/status_updates.log
    }

    # Function to mark complete
    mark_complete() {
        python3 status_utils.py complete 2>&1 | tee -a logs/status_updates.log
    }

    # Function to handle errors
    handle_error() {
        error_msg=$1
        python3 status_utils.py error "$error_msg" 2>&1 | tee -a logs/status_updates.log
        exit 1
    }

    # Monitor output for different phases
    monitor_output() {
        line_count=0
        last_update_time=$(date +%s)
        last_step=-1
        
        while IFS= read -r line; do
            echo "$line" | tee -a logs/preparation.log
            ((line_count++))

            # Convert line to lowercase for easier pattern matching
            lower_line=$(echo "$line" | tr '[:upper:]' '[:lower:]')
            
            # Detect step progression based on actual output patterns
            step_detected=0
            step_msg=""
            
            # Dynamic step detection based on configured workflow steps
            # Step indices depend on whether MEMEMBED is included
            
            # MEMEMBED/Orientation phase (only if not pre-oriented)
            if [[ "$lower_line" == *"memembed"* ]] || [[ "$lower_line" == *"orientation"* ]] || [[ "$lower_line" == *"embed"* ]]; then
                # Only count MEMEMBED if it's in our steps array
                for i in "${!STEPS[@]}"; do
                    if [ "${STEPS[i]}" = "MEMEMBED" ]; then
                        step_detected=$((i))
                        step_msg="MEMEMBED"
                        break
                    fi
                done
            # Packmol phase  
            elif [[ "$lower_line" == *"packmol"* ]] || [[ "$lower_line" == *"packing"* ]] || [[ "$lower_line" == *"pack"* ]]; then
                for i in "${!STEPS[@]}"; do
                    if [ "${STEPS[i]}" = "Packmol" ]; then
                        step_detected=$((i))
                        step_msg="Packmol"
                        break
                    fi
                done
            # pdb4amber phase
            elif [[ "$lower_line" == *"pdb4amber"* ]] || [[ "$lower_line" == *"preparing pdb"* ]]; then
                for i in "${!STEPS[@]}"; do
                    if [ "${STEPS[i]}" = "pdb4amber" ]; then
                        step_detected=$((i))
                        step_msg="pdb4amber"
                        break
                    fi
                done
            # tleap phase
            elif [[ "$lower_line" == *"tleap"* ]] || [[ "$lower_line" == *"leap"* ]] || [[ "$lower_line" == *"amber"* ]]; then
                for i in "${!STEPS[@]}"; do
                    if [ "${STEPS[i]}" = "tleap" ]; then
                        step_detected=$((i))
                        step_msg="tleap"
                        break
                    fi
                done
            fi
            
            # Update status if we detected a new step
            if (( step_detected > last_step )); then
                update_status $step_detected "$step_msg"
                last_step=$step_detected
                last_update_time=$(date +%s)
            fi
            
            # Fallback: Update based on line count and time (every 100 lines or 45 seconds)
            current_time=$(date +%s)
            elapsed=$((current_time - last_update_time))
            
            if (( line_count % 100 == 0 || elapsed > 45 )); then
                if (( last_step < 1 && line_count > 50 )); then
                    update_status 1 "Processing"
                    last_step=1
                elif (( last_step < 2 && line_count > 200 )); then
                    update_status 2 "Building"
                    last_step=2
                elif (( last_step < 3 && line_count > 400 )); then
                    update_status 3 "Finalizing"
                    last_step=3
                fi
                last_update_time=$current_time
            fi
        done
    }

    # Main execution
    echo "Starting membrane preparation at $(date)" | tee -a logs/preparation.log
    echo "Command: packmol-memgen --pdb protein_protonated_prepared.pdb --lipids POPC//POPC --ratio 1//1 --ffwat tip3p --ffprot ff14SB --fflip lipid21 --preoriented --parametrize --salt --saltcon 0.15 --dist_wat 17.5" | tee -a logs/preparation.log
    
    # Initial progress update - ensure it works
    echo "Updating initial status..." | tee -a logs/preparation.log
    update_status 0 "Starting"
    
    # Brief delay to ensure status file is written
    sleep 1

    # Execute command with monitoring and capture exit status
    echo "Launching main command..." | tee -a logs/preparation.log
    packmol-memgen --pdb protein_protonated_prepared.pdb --lipids POPC//POPC --ratio 1//1 --ffwat tip3p --ffprot ff14SB --fflip lipid21 --preoriented --parametrize --salt --saltcon 0.15 --dist_wat 17.5 2>&1 | monitor_output &
    CMD_PID=$!
    
    # Enhanced completion checker - run in background
    {
        echo "Background completion checker started" >> logs/preparation.log
        sleep 30  # Wait 30 seconds before first check
        
        while kill -0 $CMD_PID 2>/dev/null; do
            sleep 15  # Check every 15 seconds
            
            # Check if output files exist to determine completion
            output_files_found=0
            for pattern in "bilayer_*.pdb" "*.pdb" "*.prmtop" "*.inpcrd"; do
                if compgen -G "$pattern" > /dev/null 2>&1; then
                    output_files_found=1
                    echo "Found output files matching $pattern - job progressing" >> logs/preparation.log
                    break
                fi
            done
            
            if [ $output_files_found -eq 1 ]; then
                # Update to Packmol step when files are first detected (still in packing phase)
                update_status 2 "Packmol"
            fi
        done
        echo "Background completion checker finished" >> logs/preparation.log
    } &
    COMPLETION_CHECKER_PID=$!
    
    # Wait for main command to complete
    wait $CMD_PID
    EXIT_STATUS=$?

    # Clean up background checker
    if kill -0 $COMPLETION_CHECKER_PID 2>/dev/null; then
        kill $COMPLETION_CHECKER_PID 2>/dev/null || true
    fi

    # Check exit status and handle completion
    echo "Command completed with exit status $EXIT_STATUS" | tee -a logs/preparation.log
    
    if [ $EXIT_STATUS -eq 0 ]; then
        echo "Membrane preparation completed successfully at $(date)" | tee -a logs/preparation.log
        
        # Check if parametrization was requested
        parametrize_requested=true
        
        if [ "$parametrize_requested" = "true" ]; then
            echo "Starting parametrization workflow..." | tee -a logs/preparation.log
            
            # Find the bilayer PDB file generated by packmol-memgen
            bilayer_pdb=""
            for pdb_file in bilayer_*.pdb; do
                if [ -f "$pdb_file" ]; then
                    bilayer_pdb="$pdb_file"
                    echo "Found bilayer PDB file: $bilayer_pdb" | tee -a logs/preparation.log
                    break
                fi
            done
            
            if [ -z "$bilayer_pdb" ]; then
                echo "❌ No bilayer PDB file found for parametrization" | tee -a logs/preparation.log
                handle_error "Bilayer PDB file not found"
            fi
            
            # Step 1: Run pdb4amber
            echo "Running pdb4amber on $bilayer_pdb..." | tee -a logs/preparation.log
            
            # Find the correct step index for pdb4amber in the dynamic steps array
            pdb4amber_step_index=-1
            for i in "${!STEPS[@]}"; do
                if [ "${STEPS[i]}" = "pdb4amber" ]; then
                    pdb4amber_step_index=$i
                    break
                fi
            done
            
            if [ $pdb4amber_step_index -ge 0 ]; then
                update_status $pdb4amber_step_index "pdb4amber"
            else
                echo "Warning: pdb4amber step not found in steps array" | tee -a logs/preparation.log
            fi
            
            prepared_pdb="system_for_tleap.pdb"
            if pdb4amber -i "$bilayer_pdb" -o "$prepared_pdb" 2>&1 | tee -a logs/parametrization.log; then
                echo "✅ pdb4amber completed successfully" | tee -a logs/preparation.log
                
                # Step 2: Run tleap
                echo "Running tleap parametrization..." | tee -a logs/preparation.log
                
                # Find the correct step index for tleap in the dynamic steps array
                tleap_step_index=-1
                for i in "${!STEPS[@]}"; do
                    if [ "${STEPS[i]}" = "tleap" ]; then
                        tleap_step_index=$i
                        break
                    fi
                done
                
                if [ $tleap_step_index -ge 0 ]; then
                    update_status $tleap_step_index "tleap"
                else
                    echo "Warning: tleap step not found in steps array" | tee -a logs/preparation.log
                fi
                
                # Get force field configuration
                protein_ff="ff14SB"
                lipid_ff="lipid21"
                water_model="tip3p"
                
                # Map force fields to leaprc files
                case "$protein_ff" in
                    "ff14SB") protein_leaprc="leaprc.protein.ff14SB" ;;
                    "ff19SB") protein_leaprc="leaprc.protein.ff19SB" ;;
                    "ff99SB") protein_leaprc="leaprc.protein.ff99SB" ;;
                    "ff03") protein_leaprc="leaprc.protein.ff03" ;;
                    *) protein_leaprc="leaprc.protein.ff14SB" ;;
                esac
                
                case "$lipid_ff" in
                    "lipid21") lipid_leaprc="leaprc.lipid21" ;;
                    "lipid17") lipid_leaprc="leaprc.lipid17" ;;
                    "GAFF") lipid_leaprc="leaprc.gaff" ;;
                    *) lipid_leaprc="leaprc.lipid21" ;;
                esac
                
                case "$water_model" in
                    "tip3p") water_leaprc="leaprc.water.tip3p" ;;
                    "tip4p") water_leaprc="leaprc.water.tip4pew" ;;
                    "spce") water_leaprc="leaprc.water.spce" ;;
                    "opc") water_leaprc="leaprc.water.opc" ;;
                    *) water_leaprc="leaprc.water.tip3p" ;;
                esac
                
                echo "Force field configuration:" | tee -a logs/preparation.log
                echo "  Protein FF: $protein_ff -> $protein_leaprc" | tee -a logs/preparation.log
                echo "  Lipid FF: $lipid_ff -> $lipid_leaprc" | tee -a logs/preparation.log
                echo "  Water model: $water_model -> $water_leaprc" | tee -a logs/preparation.log
                
                # Create dynamic tleap input file
                cat > leap_parametrize.in << EOF
# Load force field for proteins $protein_ff
source $protein_leaprc

# Load force field for lipids $lipid_ff
source $lipid_leaprc

# Load water model ${water_model^^}
source $water_leaprc

# Load PDB file prepared by pdb4amber (protein + membrane + water + neutralized)
system = loadPDB PREPARED_PDB_PLACEHOLDER

# Check total system charge
charge system

# Neutralize total charge
addIonsRand system Na+ 0
addIonsRand system Cl- 0

# Save parameter and coordinate files
saveAmberParm system system.prmtop system.inpcrd

# Save the system processed by tleap as PDB
savePDB system system.pdb

# Exit
quit
EOF
                
                # Replace placeholder with actual file path
                sed -i "s|PREPARED_PDB_PLACEHOLDER|$(pwd)/$prepared_pdb|g" leap_parametrize.in
                
                # Execute tleap
                if tleap -f leap_parametrize.in 2>&1 | tee -a logs/parametrization.log; then
                    echo "✅ tleap parametrization completed successfully" | tee -a logs/preparation.log
                    
                    # Check if output files were created
                    output_files=("system.prmtop" "system.inpcrd" "system.pdb")
                    all_created=true
                    
                    echo "Generated files:" | tee -a logs/preparation.log
                    for filename in "${output_files[@]}"; do
                        if [ -f "$filename" ]; then
                            size=$(stat -c%s "$filename" 2>/dev/null || stat -f%z "$filename" 2>/dev/null || echo "unknown")
                            echo "   $filename: $size bytes" | tee -a logs/preparation.log
                        else
                            echo "   ❌ $filename: NOT CREATED" | tee -a logs/preparation.log
                            all_created=false
                        fi
                    done
                    
                    if [ "$all_created" = true ]; then
                        echo "✅ All parametrization files created successfully" | tee -a logs/preparation.log
                    else
                        echo "❌ Some parametrization files were not created" | tee -a logs/preparation.log
                        handle_error "Incomplete parametrization output"
                    fi
                    
                else
                    echo "❌ tleap parametrization failed" | tee -a logs/preparation.log
                    handle_error "tleap parametrization failed"
                fi
                
            else
                echo "❌ pdb4amber failed" | tee -a logs/preparation.log
                handle_error "pdb4amber failed"
            fi
            
        else
            echo "Skipping parametrization - not requested" | tee -a logs/preparation.log
        fi
        
        # Final status update
        update_status 5 "Completed"
        
        # Mark as completed 
        mark_complete
        
    else
        echo "Command failed with exit status $EXIT_STATUS" | tee -a logs/preparation.log
        
        # Check if output files exist anyway (sometimes commands exit with error but produce output)
        output_found=0
        for pattern in "bilayer_*.pdb" "*.pdb" "*.prmtop" "*.inpcrd"; do
            if compgen -G "$pattern" > /dev/null 2>&1; then
                output_found=1
                echo "Found output files despite exit status - marking as completed" | tee -a logs/preparation.log
                update_status 5 "Completed"
                mark_complete
                break
            fi
        done
        
        if [ $output_found -eq 0 ]; then
            echo "No output files found - marking as failed" | tee -a logs/preparation.log
            handle_error "Preparation failed with exit code $EXIT_STATUS"
        fi
    fi
    