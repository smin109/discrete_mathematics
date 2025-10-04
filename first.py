TOL = 1e-12

def read_matrix():
    n = int(input("n = "))
    print("행렬을 입력하세요 (행 단위, 공백 구분)")
    A = [list(map(float, input().split())) for _ in range(n)]
    return A

def I(n):
    return [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]

def det(A):
    n = len(A)
    M = [r[:] for r in A]
    d, s = 1.0, 1.0
    for c in range(n):
        p = max(range(c, n), key=lambda r: abs(M[r][c]))
        if abs(M[p][c]) < TOL: return 0.0
        if p != c:
            M[c], M[p] = M[p], M[c]
            s = -s
        pv = M[c][c]; d *= pv
        for r in range(c+1, n):
            f = M[r][c]/pv
            for k in range(c, n):
                M[r][k] -= f*M[c][k]
    return s*d

def minorM(A, r, c):
    return [row[:c]+row[c+1:] for i,row in enumerate(A) if i!=r]

def inv_adjugate(A):
    n = len(A)
    d = det(A)
    if abs(d) < TOL: return None, d
    C = [[0.0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            m = det(minorM(A, i, j))
            C[i][j] = -m if (i+j)%2 else m
    adj = [[C[j][i] for j in range(n)] for i in range(n)]
    return [[adj[i][j]/d for j in range(n)] for i in range(n)], d

def copyM(A):
    return [r[:] for r in A]

def inv_gj(A):
    n = len(A)
    M, B = copyM(A), I(n)
    for c in range(n):
        p = max(range(c, n), key=lambda r: abs(M[r][c]))
        if abs(M[p][c]) < TOL: return None
        if p != c:
            M[c], M[p] = M[p], M[c]
            B[c], B[p] = B[p], B[c]
        pv = M[c][c]
        for k in range(n):
            M[c][k] /= pv; B[c][k] /= pv
        for r in range(n):
            if r == c: continue
            f = M[r][c]
            for k in range(n):
                M[r][k] -= f*M[c][k]
                B[r][k] -= f*B[c][k]
    return B

def same(A, B, tol=1e-8):
    if A is None or B is None: return False
    n, m = len(A), len(A[0])
    for i in range(n):
        for j in range(m):
            if abs(A[i][j]-B[i][j]) > tol: return False
    return True

def show(M, title):
    print(f"\n{title}:")
    if M is None:
        print("없음")
        return
    for r in M:
        print(" ".join(f"{x:.6f}" for x in r))

def main():
    A = read_matrix()
    inv1, d = inv_adjugate(A)
    show(inv1, "[행렬식(수반행렬) 역행렬]")
    inv2 = inv_gj(A)
    show(inv2, "[가우스-조던 역행렬]")
    if inv1 is None or inv2 is None:
        print("\n[비교] 하나 이상 없음 → 비교 불가")
    else:
        print("\n[비교] 동일한가? →", "예" if same(inv1, inv2) else "아니오")
    print(f"\n참고 det(A) = {d:.6f}")

if __name__ == "__main__":
    main()
