struct Stack
{
	int mData[1000];
	int mLen;
};
//��ʼ��ջ
void InitStack(Stack &S)
{
	S.mLen = 0;
}
//Ԫ�ؽ�ջ
void Push(Stack &S,int item)
{
	S.mData[S.mLen++] = item;
}
//ɾ��ջ��Ԫ��
int Pop(Stack &S)
{
	S.mLen--;
	return S.mData[S.mLen];
}
//����ջ��Ԫ��
int  Peek(Stack &S)
{
	return S.mData[S.mLen-1];
}
//�ж�ջ�Ƿ�Ϊ��
bool EmptyStack(Stack &S)
{
	if(S.mLen == 0) return true;
	return false;
}
//���ջ
void Clear(Stack &S)
{
	for(int i = 0;i<S.mLen;++i)
	{
		Pop(S);
	}
}
