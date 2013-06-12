#include <stdio.h>
#include <string.h>
#ifdef WIN32
#include <winsock2.h>
#define NEWLINE "\r\n"
#else
#include <netinet/in.h>
#include <sys/socket.h>
#define NEWLINE "\n"
#ifndef INVALID_SOCKET
#define INVALID_SOCKET  -1
#endif
#ifndef SOCKET_ERROR
#define SOCKET_ERROR  -1
#endif
#endif


int main(int argc, char **argv)
{
  char command[1000];
  int port_start = 49152;
  int error = 0;
  struct sockaddr_in client_sockaddr;
  #ifdef WIN32
  WORD wVersionRequested;
  WSADATA wsaData;
  SOCKET sock;
  #else
  int sock;
  #endif

  
  #ifdef WIN32
  wVersionRequested = MAKEWORD(2, 2);

  error = WSAStartup(wVersionRequested, &wsaData);
  if (!error) {
    error = ((LOBYTE(wsaData.wVersion) != 2)
      || (HIBYTE(wsaData.wVersion) != 2));
    if (error) {
      WSACleanup();
    }
  }
  if (error) {
    return -1;
  }
  #endif
  error = 0;
  while ((!error) && (port_start < 65535)) {
    sock = socket(AF_INET, SOCK_STREAM, 0);
    if (sock == INVALID_SOCKET) {
      printf("INVALID_SOCKET\n");
      #ifdef WIN32
      WSACleanup();
      #endif
      return -1;
    }
    memset(&client_sockaddr, 0, sizeof(client_sockaddr));
    client_sockaddr.sin_family = AF_INET;
    client_sockaddr.sin_addr.s_addr = inet_addr(LOCALHOST_IP);
    client_sockaddr.sin_port = htons(port_start);
    printf("trying port = %d\n", port_start);
    error = connect(sock, (struct sockaddr *)&client_sockaddr,
      sizeof(client_sockaddr));
    if (!error) {
      ++port_start;
      #ifdef WIN32
      closesocket(sock);
      #else
      close(socket);
      #endif
    }
  }
  if (error) {
    printf("OK, port = %d\n", port_start);
  }
  #ifdef WIN32
  closesocket(sock);
  WSACleanup();
  #else
  close(sock);
  #endif
  
  return 0;
}
