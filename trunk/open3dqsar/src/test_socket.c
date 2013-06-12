#include <stdio.h>
#include <string.h>
#define LOCALHOST_IP "127.0.0.1"
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
  int port_num[] = {13000, 14540, 26201, 25666, 16280, 20047, 20091, 10849, 19441, 29408, 29449};
  int error = 0;
  int i = 0;
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
  sock = socket(AF_INET, SOCK_STREAM, 0);
  if (sock == INVALID_SOCKET) {
    printf("INVALID_SOCKET\n");
    #ifdef WIN32
    WSACleanup();
    #endif
    return -1;
  }
  error = SOCKET_ERROR;
  while ((error == SOCKET_ERROR) && (i < 10)) {
    client_sockaddr.sin_family = AF_INET;
    client_sockaddr.sin_addr.s_addr = inet_addr(LOCALHOST_IP);
    client_sockaddr.sin_port = htons(port_num[i]);
    printf("port = %d\n", port_num[i]);
    error = connect(sock, (struct sockaddr *)&client_sockaddr,
      sizeof(client_sockaddr));
    ++i;
  }
  if (error == SOCKET_ERROR) {
    printf("1) SOCKET_ERROR\n");
    #ifdef WIN32
    closesocket(sock);
    WSACleanup();
    #else
    close(sock);
    #endif
    return -1;
  }
  sprintf(command, "%s%s%s%s", "{\"type\":\"command\", \"command\":",
    "\"background white;background red\"", "}", NEWLINE);
  printf("command = '%s'\n", command);
  error = send(sock, command, strlen(command), 0);
  if (error == SOCKET_ERROR) {
    printf("2) SOCKET_ERROR\n");
    #ifdef WIN32
    closesocket(sock);
    WSACleanup();
    #else
    close(sock);
    #endif
    return -1;
  }
  printf("OK\n");
  #ifdef WIN32
  closesocket(sock);
  WSACleanup();
  #else
  close(sock);
  #endif
  
  return 0;
}
