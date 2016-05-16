struct Stack
{
	int mData[1000];
	int mLen;
};
//初始化栈
void InitStack(Stack &S)
{
	S.mLen = 0;
}
//元素进栈
void Push(Stack &S,int item)
{
	S.mData[S.mLen++] = item;
}
//删除栈顶元素
int Pop(Stack &S)
{
	S.mLen--;
	return S.mData[S.mLen];
}
//返回栈顶元素
int  Peek(Stack &S)
{
	return S.mData[S.mLen-1];
}
//判断栈是否为空
bool EmptyStack(Stack &S)
{
	if(S.mLen == 0) return true;
	return false;
}
//清空栈
void Clear(Stack &S)
{
	for(int i = 0;i<S.mLen;++i)
	{
		Pop(S);
	}
}
