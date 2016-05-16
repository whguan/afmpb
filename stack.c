//**************************************************************************************
// iterative with stack
//***************************************************************************************
struct stack{
    elemType *stack;        /* 存储栈元素的数组指针 */
    int top;                /* 存储栈顶元素的下标位置 */
    int maxSize;            /* 存储stack数组的长度 */
};

void againMalloc(struct stack *s)
{
    /* 空间扩展为原来的2倍，原内容被自动拷贝到p所指向的存储空间中 */
    elemType *p;
    p = realloc(s->stack, 2 * s->maxSize * sizeof(elemType));
    if(!p){
        printf("内在分配失败！");
        exit(1);
    }
    s->stack = p;        /* 使stack指向新的栈空间 */
    s->maxSize = 2 * s->maxSize;        /* 把栈空间的大小修改新的长度 */
    return;
}

/* 1.初始化栈s为空 */
void initStack(struct stack *s, int ms)
{
    if(ms <=0){
        printf("ms的值非法！");
        exit(1);
    }
    s->maxSize = ms;
    s->stack = malloc(ms * (sizeof(elemType)));
    if(!s->stack){
        printf("内在分配失败！");
        exit(1);
    }
    s->top = -1;        /* 初始置栈为空 */
    return;
}

/* 2.新元素进栈，即把它插入到栈顶 */
void push(struct stack *s, elemType x)
{
    /* 若栈空间用尽则重新分配更大的存储空间 */
    if(s->top == s->maxSize - 1){
        againMalloc(s);
    }
    s->top++;        /* 栈顶指针后移一个位置 */
    s->stack[s->top] = x;        /* 将新元素插入到栈顶 */
    return;
}

/* 3.删除(弹出)栈顶元素并返回其值 */
elemType pop(struct stack *s)
{
    /* 若栈空则退出运行 */
    if(s->top == -1){
        printf("栈空，无元素出栈！");
        exit(1);
    }
    s->top--;        /* 栈顶指针减1表示出栈 */
    return s->stack[s->top+1];        /* 返回原栈顶元素的值 */
}

/* 4.读取栈顶元素的值 */
elemType peek(struct stack *s)
{
    /* 若栈空则退出运行 */
    if(s->top == -1){
        printf("栈空，无元素可读取！");
        exit(1);
    }
    return s->stack[s->top];        /* 返回原栈顶元素的值 */
}

/* 5.判断s是否为空，若是则返回1表示真，否则返回0表示假 */
int emptyStack(struct stack *s)
{
    if(s->top == -1){
        return 1;
    }else{
        return 0;
    }
}

/* 6.清除栈s中的所有元素，释放动态存储空间 */
void clearStack(struct stack *s)
{
    if(s->stack){
        free(s->stack);
        s->stack = NULL;
        s->top = -1;
        s->maxSize = 0;
    }
    return;
}




void AggregateSweep(const int ibox)
{   // run a post-order traversal on a tree
    // use visited[x]  to record whether the children of x have been pushed onto the stask
    Stack Q;
    Q.push(ibox);
    bool visited[ibox] = {false};
    while (! Q.isEmpty() )
        t=Q.peek();
        if (g_sboxes[t].nchild == 0) {  // t is a leaf
                                                                                                       155,4         50%
        s->top = -1;
        s->maxSize = 0;
    }
    return;
}




void AggregateSweep(const int ibox)
{   // run a post-order traversal on a tree
    // use visited[x]  to record whether the children of x have been pushed onto the stask
    struct stack Q; 
    initStack(&Q,ibox)
    push(&s,ibox);
    bool visited[ibox] = {false};
    while (! Q.isEmpty() )
        t=Q.peek();
        if (g_sboxes[t].nchild == 0) {  // t is a leaf
            SourceToMultipole(t);
            if (t > 1)
            MultipoleToExponential(t);
            Q.pop(t);        } else   {                      //t has children            if visited[t]==false {
                cilk_for (int i = 0; i < 8; i++) {
                int child = g_sboxes[t].child[i];
                if (child) {
                    Q.push(child);
                    visited[child]=false;
                  }
                }            visited[t]=true;
            } else {
                MultipoleToMultipole(t);
                if (t > 1)
                MultipoleToExponential(t);
                Q.pop(t);            }
        }
}
//*************************************************************************************
// iterative//*************************************************************************************
/*struct recordStruct{        int iboxNum;//用来记录要访问节点在g_sboxes中的位置        int childNum;//用来记录下一个要访问的当前节点的子节点，不管存不存在
};
void AggregateSweep(const int ibox){
        struct recordStruct stack[1000];
                                193,5-12      64%
                    Q.push(child);
                    visited[child]=false;
                  }
                }
            visited[t]=true;
            } else {
                MultipoleToMultipole(t);
                if (t > 1)
                MultipoleToExponential(t);
                Q.pop(t);
            }
        }
}

//*************************************************************************************
// iterative
//*************************************************************************************
/*struct recordStruct{
        int iboxNum;//用来记录要访问节点在g_sboxes中的位置
        int childNum;//用来记录下一个要访问的当前节点的子节点，不管存不存在
};
void AggregateSweep(const int ibox){
        struct recordStruct stack[1000];
        //初始化
        int top = -1;
        stack[++top].iboxNum = ibox;
        stack[top].childNum = 0;
        //开始遍历
        while(top > -1 ){
                struct recordStruct temp = stack[top];//取栈顶元素
                int i = temp.childNum;
                while(i < 8){//看当前节点有没有子节点还没有访问
                        if(g_sboxes[temp.iboxNum].child[i]){//如果有，入栈
                                stack[++top].iboxNum = g_sboxes[temp.iboxNum].child[i];
                                stack[top].childNum = 0;
                                temp.childNum = ++i;//修改当前节点下一个要访问的子节点序号
                                break;
                        }//end if
                        i++;
                }//end while
                if(i == 8 ){//如果当前节点没有子节点或者子节点都访问过了，就访问当前节点并退栈
                        SourceToMultipole(temp.iboxNum);
                        if(temp.iboxNum > 1)
                                MultipoleToExponential(temp.iboxNum);
                        top--;
                        //temp = NULL;
                }
        }//end while
}
*/

