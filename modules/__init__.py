"""Helper functions for the SYBR pipeline modules"""

import os
import glob

def should_run_rule(condition_key, config):
    """Check if a rule should run based on configuration"""
    return config.get("run_stages", {}).get(condition_key, True)
