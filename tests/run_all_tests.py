import os

print("=========================== running unit tests ========================")
os.system("testflo unit_tests -v")

print("======================== running regression tests =====================")
os.system("testflo reg_tests -v")
