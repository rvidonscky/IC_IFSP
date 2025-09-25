#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Welcome Script for IFSP Scientific Initiation Repository
--------------------------------------------------------
This script displays a welcome message for visitors and contributors
to the IFSP scientific initiation projects repository.
"""

import textwrap


def main():
    message = textwrap.dedent("""
    ====================================================
               👋 Welcome to the Repository!            
    ====================================================
    
    This is the official GitHub repository for 
    IFSP Scientific Initiation Students.
    
    🚀 Here you will find:
       • Research codes
       • Experiments
       • Documentation
       • Learning resources
    
    💡 Feel free to explore, learn, and contribute!
    
    ====================================================
                Happy Coding & Research!               
    ====================================================
    """)
    print(message)


if __name__ == "__main__":
    main()
