//**************************************************************************************
// iterative with stack
//***************************************************************************************
struct stack{
    elemType *stack;        /* �洢ջԪ�ص�����ָ�� */
    int top;                /* �洢ջ��Ԫ�ص��±�λ�� */
    int maxSize;            /* �洢stack����ĳ��� */
};

void againMalloc(struct stack *s)
{
    /* �ռ���չΪԭ����2����ԭ���ݱ��Զ�������p��ָ��Ĵ洢�ռ��� */
    elemType *p;
    p = realloc(s->stack, 2 * s->maxSize * sizeof(elemType));
    if(!p){
        printf("���ڷ���ʧ�ܣ�");
        exit(1);
    }
    s->stack = p;        /* ʹstackָ���µ�ջ�ռ� */
    s->maxSize = 2 * s->maxSize;        /* ��ջ�ռ�Ĵ�С�޸��µĳ��� */
    return;
}

/* 1.��ʼ��ջsΪ�� */
void initStack(struct stack *s, int ms)
{
    if(ms <=0){
        printf("ms��ֵ�Ƿ���");
        exit(1);
    }
    s->maxSize = ms;
    s->stack = malloc(ms * (sizeof(elemType)));
    if(!s->stack){
        printf("���ڷ���ʧ�ܣ�");
        exit(1);
    }
    s->top = -1;        /* ��ʼ��ջΪ�� */
    return;
}

/* 2.��Ԫ�ؽ�ջ�����������뵽ջ�� */
void push(struct stack *s, elemType x)
{
    /* ��ջ�ռ��þ������·������Ĵ洢�ռ� */
    if(s->top == s->maxSize - 1){
        againMalloc(s);
    }
    s->top++;        /* ջ��ָ�����һ��λ�� */
    s->stack[s->top] = x;        /* ����Ԫ�ز��뵽ջ�� */
    return;
}

/* 3.ɾ��(����)ջ��Ԫ�ز�������ֵ */
elemType pop(struct stack *s)
{
    /* ��ջ�����˳����� */
    if(s->top == -1){
        printf("ջ�գ���Ԫ�س�ջ��");
        exit(1);
    }
    s->top--;        /* ջ��ָ���1��ʾ��ջ */
    return s->stack[s->top+1];        /* ����ԭջ��Ԫ�ص�ֵ */
}

/* 4.��ȡջ��Ԫ�ص�ֵ */
elemType peek(struct stack *s)
{
    /* ��ջ�����˳����� */
    if(s->top == -1){
        printf("ջ�գ���Ԫ�ؿɶ�ȡ��");
        exit(1);
    }
    return s->stack[s->top];        /* ����ԭջ��Ԫ�ص�ֵ */
}

/* 5.�ж�s�Ƿ�Ϊ�գ������򷵻�1��ʾ�棬���򷵻�0��ʾ�� */
int emptyStack(struct stack *s)
{
    if(s->top == -1){
        return 1;
    }else{
        return 0;
    }
}

/* 6.���ջs�е�����Ԫ�أ��ͷŶ�̬�洢�ռ� */
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
/*struct recordStruct{        int iboxNum;//������¼Ҫ���ʽڵ���g_sboxes�е�λ��        int childNum;//������¼��һ��Ҫ���ʵĵ�ǰ�ڵ���ӽڵ㣬���ܴ治����
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
        int iboxNum;//������¼Ҫ���ʽڵ���g_sboxes�е�λ��
        int childNum;//������¼��һ��Ҫ���ʵĵ�ǰ�ڵ���ӽڵ㣬���ܴ治����
};
void AggregateSweep(const int ibox){
        struct recordStruct stack[1000];
        //��ʼ��
        int top = -1;
        stack[++top].iboxNum = ibox;
        stack[top].childNum = 0;
        //��ʼ����
        while(top > -1 ){
                struct recordStruct temp = stack[top];//ȡջ��Ԫ��
                int i = temp.childNum;
                while(i < 8){//����ǰ�ڵ���û���ӽڵ㻹û�з���
                        if(g_sboxes[temp.iboxNum].child[i]){//����У���ջ
                                stack[++top].iboxNum = g_sboxes[temp.iboxNum].child[i];
                                stack[top].childNum = 0;
                                temp.childNum = ++i;//�޸ĵ�ǰ�ڵ���һ��Ҫ���ʵ��ӽڵ����
                                break;
                        }//end if
                        i++;
                }//end while
                if(i == 8 ){//�����ǰ�ڵ�û���ӽڵ�����ӽڵ㶼���ʹ��ˣ��ͷ��ʵ�ǰ�ڵ㲢��ջ
                        SourceToMultipole(temp.iboxNum);
                        if(temp.iboxNum > 1)
                                MultipoleToExponential(temp.iboxNum);
                        top--;
                        //temp = NULL;
                }
        }//end while
}
*/

