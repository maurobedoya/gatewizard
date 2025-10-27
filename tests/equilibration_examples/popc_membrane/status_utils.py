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
