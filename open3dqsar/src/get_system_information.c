/*

get_system_information.c

is part of

Open3DQSAR
----------

An open-source software aimed at high-throughput
chemometric analysis of molecular interaction fields

Copyright (C) 2009-2014 Paolo Tosco, Thomas Balle

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.

For further information, please contact:

Paolo Tosco, PhD
Dipartimento di Scienza e Tecnologia del Farmaco
Universita' degli Studi di Torino
Via Pietro Giuria, 9
10125 Torino (Italy)
Phone:  +39 011 670 7680
Mobile: +39 348 553 7206
Fax:    +39 011 670 7687
E-mail: paolo.tosco@unito.it

*/

#include <include/o3header.h>
#ifdef WIN32
#include <windows.h>
#include <tchar.h>

typedef void (WINAPI *PGNSI)(LPSYSTEM_INFO);

BOOL GetOSDisplayString(LPTSTR pszOS, int *page_size)
{
  OSVERSIONINFOEX osvi;
  SYSTEM_INFO si;
  PGNSI pGNSI;
  BOOL bOsVersionInfoEx;

  ZeroMemory(&si, sizeof(SYSTEM_INFO));
  ZeroMemory(&osvi, sizeof(OSVERSIONINFOEX));

  osvi.dwOSVersionInfoSize = sizeof(OSVERSIONINFOEX);

  if (!(bOsVersionInfoEx = GetVersionEx((OSVERSIONINFO *)&osvi))) {
    sprintf(pszOS, "Unsupported Microsoft OS");
    return TRUE;
  }

  // Call GetNativeSystemInfo if supported or GetSystemInfo otherwise.

  pGNSI = (PGNSI) GetProcAddress
    (GetModuleHandle(TEXT("kernel32.dll")), 
    "GetNativeSystemInfo");
  if (NULL != pGNSI)
    pGNSI(&si);
  else
    GetSystemInfo(&si);

  *page_size = (int)(si.dwAllocationGranularity);
  if (VER_PLATFORM_WIN32_NT == osvi.dwPlatformId && 
    osvi.dwMajorVersion > 4) {
    strncpy(pszOS, TEXT("Microsoft "), BUF_LEN);

    // Test for the specific product.

    if (osvi.dwMajorVersion == 6) {
      if (osvi.dwMinorVersion == 0) {
        if (osvi.wProductType == VER_NT_WORKSTATION)
          strncat(pszOS, TEXT("Windows Vista "), BUF_LEN);
        else
          strncat(pszOS, TEXT("Windows Server 2008 "), BUF_LEN);
      }

      if (osvi.dwMinorVersion == 1) {
        if (osvi.wProductType == VER_NT_WORKSTATION)
          strncat(pszOS, TEXT("Windows 7 "), BUF_LEN);
        else
          strncat(pszOS, TEXT("Windows Server 2008 R2 "), BUF_LEN);
      }
    }

    if (osvi.dwMajorVersion == 5 && osvi.dwMinorVersion == 2) {
      if (GetSystemMetrics(SM_SERVERR2))
        strncat(pszOS, TEXT("Windows Server 2003 R2, "), BUF_LEN);
      else if (osvi.wSuiteMask & VER_SUITE_STORAGE_SERVER)
        strncat(pszOS, TEXT("Windows Storage Server 2003"), BUF_LEN);
      else if (osvi.wProductType == VER_NT_WORKSTATION &&
        si.wProcessorArchitecture == PROCESSOR_ARCHITECTURE_AMD64)
          strncat(pszOS, TEXT("Windows XP Professional x64 Edition"), BUF_LEN);
      else
        strncat(pszOS, TEXT("Windows Server 2003, "), BUF_LEN);

      // Test for the server type.
      if (osvi.wProductType != VER_NT_WORKSTATION) {
        if (si.wProcessorArchitecture==PROCESSOR_ARCHITECTURE_IA64) {
          if (osvi.wSuiteMask & VER_SUITE_DATACENTER)
            strncat(pszOS, TEXT("Datacenter Edition for Itanium-based Systems"), BUF_LEN);
          else if (osvi.wSuiteMask & VER_SUITE_ENTERPRISE)
            strncat(pszOS, TEXT("Enterprise Edition for Itanium-based Systems"), BUF_LEN);
        }

        else if (si.wProcessorArchitecture==PROCESSOR_ARCHITECTURE_AMD64) {
          if (osvi.wSuiteMask & VER_SUITE_DATACENTER)
            strncat(pszOS, TEXT("Datacenter x64 Edition"), BUF_LEN);
          else if (osvi.wSuiteMask & VER_SUITE_ENTERPRISE)
            strncat(pszOS, TEXT("Enterprise x64 Edition"), BUF_LEN);
          else
            strncat(pszOS, TEXT("Standard x64 Edition"), BUF_LEN);
        }

        else {
          if (osvi.wSuiteMask & VER_SUITE_COMPUTE_SERVER)
            strncat(pszOS, TEXT("Compute Cluster Edition"), BUF_LEN);
          else if (osvi.wSuiteMask & VER_SUITE_DATACENTER)
            strncat(pszOS, TEXT("Datacenter Edition"), BUF_LEN);
          else if(osvi.wSuiteMask & VER_SUITE_ENTERPRISE)
            strncat(pszOS, TEXT("Enterprise Edition"), BUF_LEN);
          else if (osvi.wSuiteMask & VER_SUITE_BLADE)
            strncat(pszOS, TEXT("Web Edition"), BUF_LEN);
          else strncat(pszOS, TEXT("Standard Edition"), BUF_LEN);
        }
      }
    }

    if (osvi.dwMajorVersion == 5 && osvi.dwMinorVersion == 1) {
      strncat(pszOS, TEXT("Windows XP "), BUF_LEN);
      if (osvi.wSuiteMask & VER_SUITE_PERSONAL)
        strncat(pszOS, TEXT("Home Edition"), BUF_LEN);
      else
        strncat(pszOS, TEXT("Professional"), BUF_LEN);
    }

    if (osvi.dwMajorVersion == 5 && osvi.dwMinorVersion == 0) {
      strncat(pszOS, TEXT("Windows 2000 "), BUF_LEN);

      if (osvi.wProductType == VER_NT_WORKSTATION)
        strncat(pszOS, TEXT("Professional"), BUF_LEN);
      else {
        if (osvi.wSuiteMask & VER_SUITE_DATACENTER)
          strncat(pszOS, TEXT("Datacenter Server"), BUF_LEN);
        else if (osvi.wSuiteMask & VER_SUITE_ENTERPRISE)
          strncat(pszOS, TEXT("Advanced Server"), BUF_LEN);
        else
          strncat(pszOS, TEXT("Server"), BUF_LEN);
      }
    }

    // Include service pack (if any) and build number.

    if (_tcslen(osvi.szCSDVersion) > 0) {
      strncat(pszOS, TEXT(" "), BUF_LEN);
      strncat(pszOS, osvi.szCSDVersion, BUF_LEN);
    }

    TCHAR buf[80];

    snprintf(buf, 80, TEXT(" (build %d)"), (int)(osvi.dwBuildNumber));
    strncat(pszOS, buf, BUF_LEN);

    if (osvi.dwMajorVersion >= 6) {
      if (si.wProcessorArchitecture==PROCESSOR_ARCHITECTURE_AMD64)
        strncat(pszOS, TEXT(", 64-bit"), BUF_LEN);
      else if (si.wProcessorArchitecture==PROCESSOR_ARCHITECTURE_INTEL)
        strncat(pszOS, TEXT(", 32-bit"), BUF_LEN);
    }
      
    return FALSE; 
  }

  else {
    sprintf(pszOS, "Unknown Microsoft OS");
    return TRUE;
  }
}
#endif


void get_system_information(O3Data *od)
{
  #ifndef WIN32
  struct utsname host_info;
  #else
  char *nodename;
  char host_info[BUF_LEN];
  #endif


  tee_printf(od, "System information\n");
  tee_printf(od, "------------------\n");
  
  #ifndef WIN32
  memset(&host_info, 0, sizeof(struct utsname));
  if (uname(&host_info) == -1) {
    tee_printf(od, "Unavailable\n");
  }
  else {
    if (host_info.nodename) {
      tee_printf(od, "Hostname: %s\n", host_info.nodename);
    }
    if (host_info.machine) {
      tee_printf(od, "%s ", host_info.machine);
    }
    if (host_info.sysname) {
      tee_printf(od, "%s ", host_info.sysname);
    }
    if (host_info.release) {
      tee_printf(od, "%s ", host_info.release);
    }
    if (host_info.version) {
      tee_printf(od, "%s", host_info.version);
    }
    tee_printf(od, "\n");
  }
  od->mmap_pagesize = (int)sysconf(_SC_PAGESIZE);
  #else
  nodename = getenv("COMPUTERNAME");
  if (nodename) {
    tee_printf(od, "Hostname: %s\n", nodename);
  }
  GetOSDisplayString(host_info, &(od->mmap_pagesize));
  tee_printf(od, "%s\n", host_info);
  #endif
  
  od->n_proc = get_number_of_procs();
  if (od->n_proc) {
    tee_printf(od, "This machine has %d CPUs.\n\n", od->n_proc);
  }
  else {
    tee_printf(od, "This machine has an unknown number of CPUs; "
      "a single CPU will be assumed. This setting "
      "may be changed with the n_cpus parameter "
      "when appropriate\n\n");
    od->n_proc = 1;
  }
}
