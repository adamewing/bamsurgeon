import sys
import subprocess
import re

def check_java():
    p = subprocess.Popen(['java', '-version'], stderr=subprocess.PIPE)
    for line in p.stderr:
        line = line.decode()

        if line.startswith('java version') or line.startswith('openjdk version'):
            return True

    return False

def check_bwa():
    p = subprocess.Popen(['bwa'], stderr=subprocess.PIPE)
    for line in p.stderr:
        line = line.decode()

        if line.startswith('Version:'):
            major, minor, sub = line.strip().split()[1].split('.')
            sub = sub.split('-')[0]
            digit_pattern = re.compile(r'\D')
            sub = list(filter(None, digit_pattern.split(sub)))[0]
            if int(major) >= 0 and int(minor) >= 7 and int(sub) >= 12:
                return True
    return False


def check_samtools():
    p = subprocess.Popen(['samtools'], stderr=subprocess.PIPE)
    for line in p.stderr:
        line = line.decode()

        if line.startswith('Version:'):
            major, minor = line.strip().split()[1].split('.')[:2]
            minor = minor.split('-')[0]
            if int(major) >= 1 and int(minor) >= 2:
                return True
    return False


def check_wgsim():
    p = subprocess.Popen(['wgsim'], stderr=subprocess.PIPE)
    for line in p.stderr:
        line = line.decode()

        if line.startswith('Version:'):
            major, minor = line.strip().split()[1].split('.')[:2]
            minor = minor.split('-')[0]
            if int(major) >= 0 and int(minor) >= 2:
                return True
    return False


def check_velvet():
    p = subprocess.Popen(['velvetg'], stdout=subprocess.PIPE)
    for line in p.stdout: 
        line = line.decode()

        if line.startswith('Version'):
            major, minor = line.strip().split()[1].split('.')[:2]
            minor = minor.split('-')[0]
            if int(major) >= 1 and int(minor) >= 2:
                return True
    return False


def check_exonerate():
    p = subprocess.Popen(['exonerate'], stdout=subprocess.PIPE)
    for line in p.stdout:
        line = line.decode()

        if line.startswith('exonerate from exonerate'):
            major, minor = line.strip().split()[-1].split('.')[:2]
            minor = minor.split('-')[0]
            if int(major) >= 2 and int(minor) >= 2:
                return True
    return False


def check_python():
    return sys.version_info >= (3, 6)


if __name__ == '__main__':
    if not check_python(): sys.exit('Dependency problem: python >= 3.6 is required')
    if not check_bwa(): sys.exit('Dependency problem: bwa >= 0.7.12 not found')
    if not check_samtools(): sys.exit('Dependency problem: samtools >= 1.2 not found')
    if not check_wgsim(): sys.exit('Dependency problem: wgsim not found (required for addsv)')
    if not check_velvet(): sys.exit('Dependency problem: velvet >= 1.2 not found (required for addsv)')
    if not check_exonerate(): sys.exit('Dependency problem: exonerate >= 2.2 not found (required for addsv)')
    if not check_java(): sys.exit('Dependency problem: java not found')

    # TODO: Check pysam version

    print('All dependencies are satisfied')