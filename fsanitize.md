一个集成在 `gcc`、`clang` 编译器中的编译指令，可以有效测试程序中的一些诸如数组越界、未定义行为等情况。

举个例子：

```cpp
#include <bits/stdc++.h>
using namespace std;

const int maxn=2e5+5,mxr=1e5,maxm=1e7+5;
int head[maxn],nxt[maxn],to[maxn],f[maxn][25],cnt,dep[maxn],rot[maxn],ans[maxn];

struct node{int ls,rs,sum,mx;}t[maxm];

void add(int x,int y){to[++cnt]=y;nxt[cnt]=head[x];head[x]=cnt;}

void dfs(int x,int fa)
{
	f[x][0]=fa;
	dep[x]=dep[fa]+1;
	for(int i=1;i<20;++i)
        f[x][i]=f[f[x][i-1]][i-1];
	for(int i=head[x];i;i=nxt[i])
	{
		if(to[i]==fa) continue;
		dfs(to[i],x);
	}
}

int lca(int x,int y)
{
	if(dep[x]<dep[y]) swap(x,y);
	for(int i=25;i>=0;i--) if(dep[f[x][i]]>=dep[y]) x=f[x][i];
	if(x==y) return x;
	for(int i=25;i>=0;i--) if(f[x][i]!=f[y][i]) x=f[x][i],y=f[y][i];
	return f[x][0];
}

void pushup(int k)
{
	if(t[t[k].ls].mx>=t[t[k].rs].mx)
        t[k].mx=t[t[k].ls].mx,t[k].sum=t[t[k].ls].sum;
    else if(t[t[k].ls].mx<t[t[k].rs].mx)
        t[k].mx=t[t[k].rs].mx,t[k].sum=t[t[k].rs].sum;
}

void update(int &rt,int l,int r,int x,int v)
{
	if(!rt) rt=++cnt;
	if(l==r) {t[rt].mx+=v,t[rt].sum=x;return;}
	int mid=(l+r)/2;
	if(x<=mid) update(t[rt].ls,l,mid,x,v);
	else update(t[rt].rs,mid+1,r,x,v);
	pushup(rt);
}

void merge(int &p,int q,int l,int r)
{
	if(!p||!q) {p=p+q;return;}
	if(l==r) {t[p].mx+=t[q].mx;return;}
	int mid=(l+r)/2;
	merge(t[p].ls,t[q].ls,l,mid),merge(t[p].rs,t[q].rs,mid+1,r);
	pushup(p);
}

void calc(int x,int fa)
{
	for(int i=head[x];i;i=nxt[i])
	{
		if(to[i]==fa) continue;
		calc(to[i],x);
		merge(rot[x],rot[to[i]],1,mxr);
	}
	if(t[rot[x]].mx>0) ans[x]=t[rot[x]].sum;
}

int main()
{
	int n,m;cin>>n>>m;
	for(int i=1;i<n;i++)
	{
		int a,b;cin>>a>>b;
		add(a,b);add(b,a);
	} 
	cnt=0;
	dfs(1,0);
	while(m--)
	{
		int x,y,z;cin>>x>>y>>z;
		int tmp=lca(x,y);
		update(rot[x],1,mxr,z,1);
        update(rot[y],1,mxr,z,1);
        update(rot[tmp],1,mxr,z,-1);
        update(rot[f[tmp][0]],1,mxr,z,-1);	
    }
    calc(1,0);
	for(int i=1;i<=n;i++) cout<<ans[i]<<endl;
	return 0;
}
```

来自[雨天的尾巴 /【模板】线段树合并](https://www.luogu.com.cn/problem/P4556)。

这段程序存在倍增数组越界，我们用 `g++ -o a a.cpp -fsanitize=address`。

测试结果：`runtime error: index 25 out of bounds for type 'int [25]'`。

剩下的之后再写。