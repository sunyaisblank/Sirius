// CRWN001A.h - Window Management
#ifndef CRWN001A_H
#define CRWN001A_H

struct GLFWwindow;

class Window {
public:
    Window(int width, int height, const char* title);
    ~Window();

    bool shouldClose() const;
    void pollEvents();
    void swapBuffers();
    GLFWwindow* getGLFWwindow() { return m_Window; }

private:
    GLFWwindow* m_Window;
};

#endif // CRWN001A_H
